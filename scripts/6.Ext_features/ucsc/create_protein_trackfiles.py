#!/usr/bin/python

"""
Combine all track information into per-protein track files, where the expectation, variance, and covariance
between tracks has been precomputed based on the analytical framework underlying PertInInt

Contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import time
import signal
import gzip
import argparse
import numpy as np
from subprocess import call
from scipy.stats import shapiro, pearsonr

data_path = '/media/vineetb/t5-vineetb/'

####################################################################################################
# PATHS TO TRACK FILES
####################################################################################################

# tuples corresponding to track definitions; first element describes subset of protein positions that
#  the track spans, second element describes the per-match-state functional weights (0 -> 1)
standard_tracks = [(data_path+'interacdome/interacdome0.3-pfam31_domsbyprot-GRCh38.txt.gz',
                    data_path+'interacdome/interacdome0.3-pfam31_domainweights-GRCh38.txt.gz'),
                   (data_path+'canbind/canbind-biolip-to-ensembl_domsbyprot-GRCh38.txt.gz',
                    data_path+'canbind/canbind-biolip-to-ensembl_domainweights-GRCh38.txt.gz'),
                   (data_path+'ucscgb/hg38alignment/exons/100way-jsdconservation_domsbyprot-GRCh38.txt.gz',
                    data_path+'ucscgb/hg38alignment/exons/100way-jsdconservation_domainweights-GRCh38.txt.gz')]

# locations of "functional" regions (that should have a positive functional weight versus regions outside)
binary_tracks = [data_path+'interacdome/pfam31-domains_domsbyprot-GRCh38.txt.gz']

# per-gene likelihoods of harboring a mutation; tracks are actually one-hot encodings of genes
whole_gene_tracks = [data_path+'diffmut/diffmut-naturalvariation_wholegene-GRCh38.txt.gz']

# per-protein-position likelihood of missense mutation (to be used when precomputing expectations/variances)
mutation_biases = data_path+'gdc/somatic_mutations/Aggregate/TCGA.Aggregate.muse.mutational_biases.tsv'


####################################################################################################
# CONSTANTS
####################################################################################################

BUILD = 'GRCh38'

# these repeat domains are often false positive hits in very long proteins and are not known to be
# involved in cancer; if a protein has 5+ of these domains, we store the aggregate track only
repetitive_domains = {
    # epidermal growth factor domain (often found in extracellular domains of membrane-bound proteins)
    'PF00008_EGF', 'PF07645_EGF_CA', 'PF12661_hEGF', 'PF12662_cEGF', 'PF06247_Plasmod_Pvs28',
    'PF14670_FXa_inhibition',  # coagulation factor highly similar to EGF domain
    'PF00041_fn3',  # fibronectin also found in extracellular proteins
    # immunoglobulin domains (found in antibodies, titin, receptor tyrosine kinases)
    'PF00047_ig', 'PF13895_Ig_2', 'PF13927_Ig_3', 'PF07679_I-set', 'PF07686_V-set',
    # spectrin is usually involved in cytoskeletal structure
    'PF00435_Spectrin',
    # nebulin is found mostly in skeletal muscle and occurs repeated in large (600-900 kDa) proteins
    'PF00880_Nebulin',
    # found in huge human-specific repeats in neurons specifically; not known to be involved in cancer
    'PF06758_DUF1220'
}

# these extremely long genes take hours (and hours) to process and are not known to be involved
# in cancer (TTN, OBSCN, NEB)
long_proteins = {'ENSG00000155657', 'ENSG00000154358', 'ENSG00000183091'}


####################################################################################################
# PROCESS POSITIONAL MUTATION BIASES
####################################################################################################

def mutational_likelihoods_overall(likelihood_file, seq_id_subset=None):
    """
    :param likelihood_file: full path to tab-delimited file of background natural variation frequencies
                        with columns sequence ID, protein length, total alleles, and protein position: AC/AN; ...
    :param seq_id_subset: sequence ID to return results for (to speed up preprocessing for testing purposes)
    :return: dictionary of sequence ID -> (protein length, tuple(positional allele frequencies))
    """

    mutation_rates = {}  # prot_id -> (int(prot_len), tuple( (int(prot_pos), float(AF)),...)

    if os.path.isfile(likelihood_file):
        variant_handle = gzip.open(likelihood_file) if likelihood_file.endswith('gz') else open(likelihood_file)
        for var_line in variant_handle:
            if var_line.startswith('#'):
                continue
            protein_id, prot_length, allele_num, prot_position_count = var_line[:-1].split('\t')[:4]
            protein_id = protein_id.split('.')[0]  # remove version
            current_rates = []  # keep track of all minor allele rates at each position

            if seq_id_subset and protein_id not in seq_id_subset:
                continue

            for position_freq in prot_position_count.split(','):  # position:AC/AN;AC/AN etc. OR position:AC
                prot_pos = int(position_freq.split(':')[0])

                rates = []
                for allele_count in position_freq.split(':')[1].split(';'):  # 1/103242;2/103244 OR 2
                    if '/' in allele_count:
                        ac = float(allele_count.split('/')[0])
                        an = float(allele_count.split('/')[1])
                    else:
                        ac = float(allele_count)
                        an = float(allele_num)
                    rates.append(ac / max(an, 1.))
                current_rates.append((prot_pos, sum(rates)))
            mutation_rates[protein_id] = (int(prot_length), tuple(current_rates))

    return mutation_rates


########################################################################################################

def mutational_likelihoods_by_protein(mutation_rates, protein_id, seq_file, bin_values=False):
    """
    :param mutation_rates: dictionary of sequence ID -> (protein length, tuple( (protein position, allele frequency)...)
    :param protein_id: sequence ID to obtain a vector for
    :param seq_file: full path to a fasta formatted file
    :param bin_values: whether to bin the allele frequencies / mutation likelihoods or use the raw values
    :return: a vector of protein length containing values to be normalized (i.e., sum to 1) indicating mutation
             likelihood at each protein position
    """

    if protein_id in mutation_rates:
        protein_length = mutation_rates[protein_id][0]

        background_mutation = [1 if bin_values else -1] * protein_length  # length of protein

        for prot_pos, af in mutation_rates[protein_id][1]:  # 0-index position and allele frequency
            if prot_pos >= protein_length:
                continue

            if bin_values:
                background_mutation[prot_pos] = (1 if af <= .001 else
                                                 (2 if af <= .005 else
                                                  (3 if af <= .01 else
                                                   (4 if af <= .25 else 5))))
            else:
                background_mutation[prot_pos] = af

        # reset the unknown positions to be the average:
        if not bin_values and -1 in background_mutation:
            average_value = [mut_rate for mut_rate in background_mutation if mut_rate != -1]
            average_value = sum(average_value) / float(len(average_value))

            for i in xrange(len(background_mutation)):
                if background_mutation[i] == -1:
                    background_mutation[i] = average_value

    else:
        protein_length = sequence_lens(seq_file, [protein_id]).get(protein_id, 1)
        background_mutation = [1] * protein_length

    return background_mutation


####################################################################################################
# PROCESS PROTEINS
####################################################################################################

def reformat_time(run_time):
    """
    :param run_time: total time elapsed in seconds
    :return: string corresponding to a properly formatted (days, hours, minutes, seconds) time
    """

    minutes, seconds = divmod(run_time, 60)
    hours, minutes = divmod(minutes, 60)
    days, hours = divmod(hours, 24)
    return ':'.join(map(lambda x: str(int(x)).zfill(2), [days, hours, minutes, seconds]))


####################################################################################################

def handler(signum, frame, print_errors=False):
    """
    :param signum: signal number passed by signal.signal
    :param frame: what to do when we reach this time
    :param print_errors: boolean indicating whether to print the signum and frame (True) or not
    :return: raise an exception
    """
    if print_errors:
        print(signum)
        print(frame)
    raise Exception("Timed out")


####################################################################################################

def sequence_lens(seq_file, seq_id_subset=None):
    """
    :param seq_file: full path to a fasta formatted file
    :param seq_id_subset: set of sequence IDs to return results for (to speed up preprocessing for testing purposes)
    :return: dictionary of sequence ID -> length of corresponding sequence
    """
    seq_lens = {}  # sequence ID -> sequence length
    seq_id = ''

    fasta_handle = gzip.open(seq_file) if seq_file.endswith('.gz') else open(seq_file)
    for fasta_line in fasta_handle:
        if fasta_line.startswith('>'):
            seq_id = fasta_line[1:-1].split()[0].split('.')[0]
            if not seq_id_subset or seq_id in seq_id_subset:
                seq_lens[seq_id] = 0

        elif not seq_id_subset or seq_id in seq_id_subset:
            seq_lens[seq_id] += len(fasta_line.strip())
    fasta_handle.close()

    return seq_lens


####################################################################################################

def find_modelable_prots(track_files, seq_file, results_file):
    """
    :param track_files: list of files where the first element on non-comment lines is an Ensembl protein ID
                        that can be modeled using one of our approaches
    :param seq_file: full path to a fasta-formatted file (in case track_files is empty) where sequence IDs
                     are Ensembl protein IDs
    :param results_file: full path to a file to write a list of modelable proteins out to
    :return: results_file upon successful write
    """

    # (1) get the protein IDs from the track_files:
    protein_ids = set()
    for intrack_file in track_files:
        if not os.path.isfile(intrack_file):
            sys.stderr.write('Could not open ' + intrack_file + '. Continuing...\n')
            continue
        track_handle = gzip.open(intrack_file, 'rt') if intrack_file.endswith('gz') else open(intrack_file)
        for track_line in track_handle:
            if track_line.startswith('#'):
                continue
            protein_ids.add(track_line.strip().split()[0].split('.')[0])
        track_handle.close()

    # (2) get the locations of all proteins:
    seq_locs = {}  # protein ID -> chromosome ID / gene ID

    fasta_handle = gzip.open(seq_file) if seq_file.endswith('gz') else open(seq_file)
    for fasta_line in fasta_handle:
        if fasta_line.startswith('>'):
            seq_id = fasta_line[1:-1].split()[0].split('.')[0]
            seq_locs[seq_id] = (
                fasta_line[fasta_line.find(BUILD + ':') + len(BUILD) + 1:].split(':')[0],
                fasta_line[fasta_line.find('gene:') + 5:].split()[0].split('.')[0]
            )
    fasta_handle.close()

    # (3) finally, write out the protein identifiers
    out_handle = gzip.open(results_file, 'wt') if results_file.endswith('gz') else open(results_file, 'w')
    for seq_id in sorted(list(protein_ids)):
        if seq_id in seq_locs:
            out_handle.write(seq_locs[seq_id][0] + '/' + seq_locs[seq_id][1] + '/' + seq_id + '\n')
    out_handle.close()

    return results_file


####################################################################################################

def read_modelable_prots(results_file):
    """
    :param results_file: full path to a file containing an Ensembl protein ID per line that can be
                         modeled by 1 or more tracks
    :return: ordered list of proteins from file
    """

    if not os.path.isfile(results_file):
        sys.stderr.write('Could not find modelable protein list in '+results_file+'!\n' +
                         'Please run '+sys.argv[0]+' --protein_list '+results_file+'\n')
        sys.exit(1)

    # read in the proteins, line by line
    results_handle = gzip.open(results_file, 'rt') if results_file.endswith('gz') else open(results_file)
    modelable_prots = set([prot_line.strip() for prot_line in results_handle])
    results_handle.close()

    return sorted(list(modelable_prots))


####################################################################################################

def find_failed_track_weights(model_prot_file, path_to_track_weights, log_file, restrict_chroms=None,
                              ignore_genes=None):
    """
    :param model_prot_file: full path to a file created with find_modelable_prots
    :param path_to_track_weights: full path to a directory with subdirectories corresponding to
                                  genes with per-protein processed track weight files within those subdirectories
    :param log_file: full path to a file to write subset of protein IDs to that should have produced
                     a track weight file but failed
    :param restrict_chroms: set of chromosomes to restrict to: [str(c) for c in range(1, 23)] + ['X', 'Y', 'MT']
    :param ignore_genes: set of genes to be ignored
    :return: error_file upon successful write
    """

    # (1) search for corresponding track weight files for each modelable protein:
    redo_track_files = set()  # set of protein IDs to rerun track file processing on

    modelable_handle = gzip.open(model_prot_file, 'rt') if model_prot_file.endswith('gz') else open(model_prot_file)
    for protein_id in modelable_handle:
        track_weight_file = path_to_track_weights + protein_id.strip() + '.trackweights.tsv'

        # skip genes on chromosomes that we are not processing at the moment
        if restrict_chroms and track_weight_file.split('/')[-3] not in restrict_chroms:
            continue

        # skip genes that we are supposed to ignore:
        if ignore_genes and True in ['/'+gn+'/' in track_weight_file for gn in ignore_genes]:
            continue

        if not os.path.isfile(track_weight_file):
            redo_track_files.add(protein_id.strip())
    modelable_handle.close()

    # (2) write out to error file
    log_handle = gzip.open(log_file, 'wt') if log_file.endswith('gz') else open(log_file, 'w')
    for protein_id in sorted(list(redo_track_files)):
        log_handle.write(protein_id+'\n')
    log_handle.close()

    return log_file


####################################################################################################
# PROCESS TRACK INPUT
####################################################################################################

def process_standard_tracks(track_files):
    """
    :param track_files: list of tuples of file names, where the first element is a tab-delimited file with
                        protein ID, domain name, matchstate -> 0-index mapping and the second element is a
                        tab-delimited file with domain name, domain qualifier, matchstate, functional score
    :return: two dictionaries: (1) protein ID -> domain name -> [{0-index: match state}], and
             (2) domain name -> domain qualifier -> {match state: functional score}
    """

    region_locs = {}  # protein ID -> domain name -> {0-index: match state}
    region_wts = {}  # domain name -> (domain qualifier, track name) -> {match state: functional score}

    for region_loc_file, region_wt_file in track_files:

        # (1) get the region locations
        loc_handle = gzip.open(region_loc_file) if region_loc_file.endswith('gz') else open(region_loc_file)
        for loc_line in loc_handle:
            if loc_line.startswith('#'):
                continue
            seq_id, domain_name, mstate_index_aa = loc_line[:-1].split('\t')[:3]
            seq_id = seq_id.split('.')[0]

            if seq_id not in region_locs:
                region_locs[seq_id] = {}
            if domain_name not in region_locs[seq_id]:
                region_locs[seq_id][domain_name] = []

            # create a new instance of this domain name in this protein
            new_instance = {}
            for mstate, index, _ in [val.split(':') for val in mstate_index_aa.split(',')]:
                new_instance[int(index)] = mstate
            region_locs[seq_id][domain_name].append(new_instance)

        loc_handle.close()

        # (2) get the corresponding functional weights
        wt_handle = gzip.open(region_wt_file) if region_wt_file.endswith('gz') else open(region_wt_file)
        for wt_line in wt_handle:
            if wt_line.startswith('#'):
                continue
            domain_name, domain_qualifier, mstate, score = wt_line[:-1].split('\t')[:4]

            if domain_name not in region_wts:
                region_wts[domain_name] = {}

            if domain_qualifier not in region_wts[domain_name]:
                region_wts[domain_name][domain_qualifier] = {}

            region_wts[domain_name][domain_qualifier][mstate] = float(score)
        wt_handle.close()

    return region_locs, region_wts


####################################################################################################

def process_binary_tracks(track_files):
    """
    :param track_files: list of file names of tab-delimited files with protein ID, domain name, track range, and
                        domain range(s)
    :return: dictionary of protein ID -> domain name -> [track_start, track_end, set({positions}, {positions}, ...)]
    """

    region_locs = {}  # protein ID -> domain name -> [track_start, track_end, set({positions}, {positions}, ...)]

    for region_loc_file in track_files:

        loc_handle = gzip.open(region_loc_file) if region_loc_file.endswith('gz') else open(region_loc_file)
        for loc_line in loc_handle:
            if loc_line.startswith('#'):
                continue
            seq_id, domain_name, track_range, domain_range = loc_line[:-1].split('\t')[:4]
            seq_id = seq_id.split('.')[0]

            track_start, track_end = map(int, track_range.split('-'))

            if seq_id not in region_locs:
                region_locs[seq_id] = {}
            if domain_name not in region_locs[seq_id]:
                region_locs[seq_id][domain_name] = [track_start, track_end, set()]

            # get the set of positions that should be positively weighted in this track
            domain_coverage = set()
            for interval in domain_range.split(','):
                if '-' in interval:
                    interval_start, interval_end = map(int, interval.split('-'))
                    domain_coverage = domain_coverage.union(set(range(interval_start, interval_end+1)))
                else:
                    domain_coverage.add(int(interval))

            # restrict to those positions that actually fall into this track (incase input was malformed):
            domain_coverage = tuple(sorted(list(domain_coverage.intersection(set(range(track_start, track_end + 1))))))

            region_locs[seq_id][domain_name][2].add(domain_coverage)

        loc_handle.close()

    return region_locs


####################################################################################################

def process_whole_gene_tracks(track_files):
    """
    :param track_files: list of tab-delimited files with protein ID, gene ID, total genes, mutational likelihood,
                        description
    :return: total genes evaluated, dictionary of protein ID -> relative mutational likelihood
    """

    gene_wts = {}
    total_genes = 0

    for gene_score_file in track_files:

        gene_handle = gzip.open(gene_score_file) if gene_score_file.endswith('gz') else open(gene_score_file)
        for gene_line in gene_handle:
            if gene_line.startswith('#'):
                continue
            seq_id, _, gene_total, mut_likelihood = gene_line[:-1].split('\t')[:4]
            seq_id = seq_id.split('.')[0]

            if int(gene_total) > total_genes:
                total_genes = int(gene_total)

            gene_wts[seq_id] = float(mut_likelihood)
        gene_handle.close()

    return total_genes, gene_wts


####################################################################################################

def process_track_name(domain_name, domain_qualifier, domain_range, aggregate=False):
    """
    :param domain_name: name of the functional region we are looking at
    :param domain_qualifier: "qualifier" for the domain we are looking at
    :param domain_range: set of 0-indexed positions that this functional region spans
    :param aggregate: does this correspond to an "aggregate" track?
    :return: the name of this particular track
    """

    if domain_name.endswith('_100way_JSD'):  # sequence conservation by JSD scores
        return 'JSD_conservation'

    elif domain_name.endswith('_ExAC_JSD'):  # minor allele frequencies
        return 'ExAC_allelefreq'

    elif not domain_name.startswith('PF') and not domain_name.endswith('_JSD'):  # homology-inferred site
        return ','.join(['Homology_binding:' + str(min(domain_range)) + '-' + str(max(domain_range)) + ':' +
                         qualifier for qualifier in domain_qualifier.split(',')])

    elif domain_name.startswith('PF') and not domain_qualifier:  # Pfam domain
        if not aggregate:
            dom_span = ','.join([str(interval_start) + '-' + str(interval_end)
                                 for interval_start, interval_end in positions_to_tuple_ranges(domain_range)])
            return domain_name + ':' + dom_span + ':complete'
        else:
            return domain_name + ':aggregate:complete'

    elif domain_name.startswith('PF') and domain_qualifier:  # InteracDome domain
        if not aggregate:
            dom_span = ','.join([str(interval_start) + '-' + str(interval_end)
                                 for interval_start, interval_end in positions_to_tuple_ranges(domain_range)])
            return ','.join([domain_name + ':' + dom_span + ':' + qualifier
                             for qualifier in domain_qualifier.split(',')])
        else:
            return ','.join([domain_name + ':aggregate:' + qualifier for qualifier in domain_qualifier.split(',')])

    else:
        if not aggregate:
            dom_span = ','.join([str(interval_start)+'-'+str(interval_end)
                                 for interval_start, interval_end in positions_to_tuple_ranges(domain_range)])
            return ','.join([domain_name + ':' + dom_span + ':' + qualifier
                             for qualifier in domain_qualifier.split(',')])
        else:
            return ','.join([domain_name + ':aggregate:' + qualifier for qualifier in domain_qualifier.split(',')])


####################################################################################################
# STEP 1: Create track weight files per protein
####################################################################################################

def add_new_weight_vector(this_vector, this_name, this_position,
                          weight_vectors, vector_names, vector_positions):
    """
    :param this_vector: vector of length L (protein length) corresponding to a new set of binding weights
    :param this_name: name of vector (e.g., PF00096_zf-C2H2_DNABASE__23-46)
    :param this_position: [start/end] position of where we will be shuffling domains over
    :param weight_vectors: existing list of vectors
    :param vector_names: corresponding existing list of vector names
    :param vector_positions: corresponding existing list of where to shuffle mutations
    :return: add the new vector to the existing lists and return IFF it does not exactly match or nearly exactly
             match any vector that is already being stored. In the case of highly redundant domains, AND identical
             binding potential weights (AGG_, NUCACID_, DNA_, DNABASE_), this is quite possible.
    """

    if len(set(this_vector)) < 2:  # No variance! Impossible to calculate pearson correlation.
        return

    for i, other_vector in enumerate(weight_vectors):

        if pearsonr(this_vector, other_vector)[0] > 0.9:

            # update the current average vector
            number_merged = vector_names.count(',') + 1  # how many existing vectors have been merged in?
            weight_vectors[i] = [((w * number_merged) + this_vector[j]) / (number_merged + 1.) for j, w in
                                 enumerate(other_vector)]

            vector_names[i] += ',' + this_name

            old_positions = vector_positions[i]

            for startn, endn in this_position:  # for each NEW interval to be merged in...
                new_positions = set()
                added_new = False

                for interval_start, interval_end in old_positions:  # for each existing interval to check against
                    if startn - 5 < interval_start < startn + 5 and endn - 5 < interval_end < endn + 5:
                        new_positions.add((min(startn, interval_start), max(endn, interval_end)))
                        added_new = True
                    else:
                        new_positions.add((interval_start, interval_end))
                if not added_new:
                    new_positions.add((startn, endn))

                old_positions = new_positions
            vector_positions[i] = old_positions
            break
    else:
        weight_vectors.append(this_vector)
        vector_names.append(this_name)
        vector_positions.append(set(this_position))


####################################################################################################

def positions_to_tuple_ranges(positions):
    """
    :param positions: list of ints representing 0-indices in a sequence
    :return: ordered list of tuples describing the start/end positions of multiple ranges
    """

    indices = sorted(map(int, list(positions)))
    ranges = []

    range_start = indices[0]  # start of current range
    range_end = indices[0]  # end of current range
    for index in indices[1:]:
        if index > range_end + 1:
            ranges.append((range_start, range_end))
            range_start = index
            range_end = index
        else:
            range_end += 1
    ranges.append((range_start, range_end))

    return ranges  # e.g., [(1,2), (5,5), (20,22)...]


####################################################################################################

def create_header(protein_id, chrom_id, gene, total_tracks, whole_gene_info, track_type1, track_type2, track_type3):
    """
    :param protein_id: sequence ID for this output file
    :param chrom_id: chromosome where protein is located
    :param gene: gene with this protein as an isoform
    :param total_tracks: total tracks found in this protein
    :param whole_gene_info: information about relative mutation rate across all genes
    :param track_type1: list of all "standard" track type files
    :param track_type2: list of all "binary" track type files
    :param track_type3: list of all "whole gene" track type files
    :return: string corresponding to the header to write to output files
    """

    header = ['# Per-protein-position functional scores for ' + str(total_tracks) + ' total tracks in ' +
              protein_id + ' (chromosome ' + chrom_id + ', gene ' + gene + ', build ' + BUILD + ')',
              '# Functional protein subregions and corresponding functional scores listed in: '] + \
             ['# ('+str(i+1)+') ' + region_loc + '\n' + '#     ' + region_wts
              for i, (region_loc, region_wts) in enumerate(track_type1)] + \
             ['# Functional protein subregions with constant scores (i.e., binary) listed in: '] + \
             ['# ('+str(i+1)+') ' + region_loc for i, region_loc in enumerate(track_type2)] + \
             ['# Whole gene mutation likelihoods listed in: '] + \
             ['# ('+str(i+1)+') ' + region_loc for i, region_loc in enumerate(track_type3)] + \
             (['# Relative Mutability & Total Genes Evaluated = ' + str(whole_gene_info[0]) + ' ' +
              str(whole_gene_info[1])] if whole_gene_info else []) + \
             ['# Background per-position mutation likelihoods listed in:',
              '# (1) ' + mutation_biases,
              '# NOTE: The expectation, variance, and covariance (between tracks) assumes ONE mutation, ',
              '#       but because of the linearity of expectations, true values are obtained by multiplying ',
              '#       these values by the total number of mutations within the enrichment interval(s)',
              '\t'.join(['protein_id', 'track_id', 'track_name', '0-index-enrichment-intervals',
                         '0-index-positive-functional-scores', 'expectation_yi', 'variance_yi', 'covariance',
                         'minimum_mutation_count', 'empirical_runtime'])]
    return '\n'.join(header)+'\n'


####################################################################################################

def process_domain_weights(test_prot_id, test_gene_id, test_chrom_id, test_prot_len, track_out_file,
                           domain_locs, domain_scores, binary_domain_locs, total_genes, whole_gene_weights):
    """
    :param test_prot_id: sequence ID to get results for
    :param test_gene_id: corresponding gene ID for the protein ID
    :param test_chrom_id: corresponding chromosome that test_prot_id is location on
    :param test_prot_len: length of sequence ID that we are processing
    :param track_out_file: full path to write track results to
    :param domain_locs: dictionary of domain name -> [{0-index: match state}, ...]
    :param domain_scores: dictionary of domain name -> (ligand type, track name) -> match state -> binding score
    :param binary_domain_locs: domain name -> [track_start, track_end, set({positions}, {positions}, ...)]
    :param total_genes: total genes with relative mutational likelihood scores
    :param whole_gene_weights: mutational likelihood score of this particular protein
    :return: remove redundancies between domain weights over proteins
    """

    # Keep track of a running list of tracks (i.e., their positions, scores, and corresponding name)
    # NOTE #1: we do not store multiple tracks with correlations of 1 to any existing functional weight track.
    # NOTE #2: these three variables will be edited by the add_new_weight_vector function and NOT explicitly here
    weight_vectors = []  # list of lists corresponding to different tracks; each track is length L (protein length)
    vector_names = []  # names of the tracks (e.g., 'JSD_conservation', 'ExAC_allelefreq'
    vector_positions = []  # positions that this track spans (usually 0->len(L)-1 but not for subregions)

    # ------------------------------------------------------------
    # STEP 1: PROCESS STANDARD TRACKS
    # ------------------------------------------------------------
    if domain_locs:
        for domain_type in domain_locs.keys():

            for domain_qualifier, matchstate_mapping in domain_scores[domain_type].items():

                total_instances = 0
                aggregate_positions = set()  # keep track of 0-index positions that this track will span (all instances)
                aggregate_vector = [0.] * test_prot_len  # mapping from 0-index position -> functional score

                # for each domain instance:
                for domain_instance in domain_locs[domain_type]:

                    # (1) get the 0-index positions that this domain instance spans:
                    this_position = [index for index in domain_instance.keys() if index in range(test_prot_len)]
                    if len(this_position) < 1:
                        continue

                    # (2) get the functional scores associated with each of these track positions
                    this_vector = [0.] * test_prot_len
                    for index, mstate in domain_instance.items():
                        if index < test_prot_len:
                            this_vector[index] = matchstate_mapping.get(mstate, 0.)

                    # (3) get the track name
                    this_name = process_track_name(domain_type, domain_qualifier, this_position)

                    # (4) add new weight vector
                    if len(set(this_vector)) < 2:
                        sys.stderr.write('No nonzero weights for ' + test_prot_id + ' on ' + this_name + '...\n')
                        continue

                    # only add this weight vector if it is NOT a member of a problematic domain family OR at the very
                    #   least, if there are fewer instances of that domain family (i.e., more likely to be
                    #   biologically relevant)
                    elif domain_type not in repetitive_domains or len(domain_locs[domain_type]) < 5:
                        add_new_weight_vector(this_vector, this_name, positions_to_tuple_ranges(this_position),
                                              weight_vectors, vector_names, vector_positions)

                    # and keep track of the aggregate mapping as well, if necessary
                    total_instances += 1
                    aggregate_positions = aggregate_positions.union(set(this_position))
                    for index, score in enumerate(this_vector):
                        aggregate_vector[index] = max(aggregate_vector[index], score)

                # And finally, for all domains in aggregate!
                if total_instances > 1 and len(aggregate_positions) > 0:
                    aggregate_name = process_track_name(domain_type, domain_qualifier, aggregate_positions, True)
                    add_new_weight_vector(aggregate_vector, aggregate_name,
                                          positions_to_tuple_ranges(aggregate_positions),
                                          weight_vectors, vector_names, vector_positions)

    # ------------------------------------------------------------
    # STEP 2: PROCESS BINARY TRACKS
    # ------------------------------------------------------------
    if binary_domain_locs:
        for domain_type in binary_domain_locs.keys():

            # keep track of "protein length" (i.e., where this binary track spans)
            track_start, track_end = binary_domain_locs[domain_type][:2]
            this_position = set([index for index in range(track_start, track_end + 1) if index in range(test_prot_len)])

            # keep track of positively-scored 0-indices in case there are multiple instances
            total_instances = 0
            aggregate_positions = set()

            for domain_positions in binary_domain_locs[domain_type][2]:

                if len(domain_positions) < 1:
                    continue

                # (1) get the functional scores associated with each track position
                this_vector = [0.] * test_prot_len
                for index in domain_positions:
                    if index < test_prot_len:
                        this_vector[index] = 0.05

                # (2) and get the name for this region
                this_name = process_track_name(domain_type, None, domain_positions)

                # (3) add new track
                if len(set(this_vector)) < 2:
                    sys.stderr.write('No nonzero weights for ' + test_prot_id + ' on ' + this_name + '...\n')
                    continue
                elif domain_type not in repetitive_domains or len(binary_domain_locs[domain_type][2]) < 5:
                    add_new_weight_vector(this_vector, this_name, positions_to_tuple_ranges(this_position),
                                          weight_vectors, vector_names, vector_positions)

                # and keep track of the aggregate mapping as well, if necessary
                total_instances += 1
                aggregate_positions = aggregate_positions.union(set(domain_positions))

            # for all the domain regions in aggregate:
            if total_instances > 1 and len(aggregate_positions) > 0:
                aggregate_name = process_track_name(domain_type, None, aggregate_positions, True)
                aggregate_vector = [0.] * test_prot_len
                for index in aggregate_positions:
                    aggregate_vector[index] = 0.05
                add_new_weight_vector(aggregate_vector, aggregate_name, positions_to_tuple_ranges(this_position),
                                      weight_vectors, vector_names, vector_positions)

    # ------------------------------------------------------------
    # STEP 3: PROCESS WHOLE GENE TRACKS & WRITE OUT TO FILE
    # ------------------------------------------------------------
    for subdir in ['/'.join(track_out_file.split('/')[:i]) for i in xrange(2, track_out_file.count('/')+1)]:
        if not os.path.isdir(subdir):
            call(['mkdir', subdir])

    out_handle = open(track_out_file, 'w')
    out_handle.write(create_header(prot_id,  # protein ID
                                   test_chrom_id,  # chromosome
                                   test_gene_id,  # gene ID
                                   len(weight_vectors),  # total tracks in this protein
                                   (whole_gene_weights, total_genes) if whole_gene_weights else None,
                                   standard_tracks, binary_tracks, whole_gene_tracks))

    for i in xrange(len(weight_vectors)):
        out_handle.write('\t'.join([prot_id, str(i + 1), vector_names[i],
                                    ','.join(map(lambda x: str(x[0]) + '-' + str(x[1]),
                                                 sorted(list(vector_positions[i])))),
                                    ','.join([str(j) + ':' + str(w) for j, w in enumerate(weight_vectors[i])
                                              if w > 0])]) + '\n')

    out_handle.close()


####################################################################################################
# STEP 2: Include expectation, variance, and covariance in preprocessed track weight files
####################################################################################################

def update_domain_weights(original_file, complete_file, mutation_rates=None):
    """
    :param original_file: full path to a tab-delimited file with (at least) columns [0] protein ID, [1] track ID,
                          [2] track name, [3] 0-index enrichment interval(s), [4] 0-index non-zero binding weights
    :param complete_file: full path to a new tab-delimited file with additional columns [5] expected value assuming
                          1 mutation, [6] variance and [7] covariances with all other tracks
    :param mutation_rates: list of protein length with values corresponding to likelihood of mutation
    :return: compute the expectation, variance, and covariances between all tracks and write out to new file
    """

    if not mutation_rates:  # index -> likelihood of mutation
        mutation_rates = {}  # no entries

    intervals = {}  # track_id -> (start,end), (start,end),...
    nonzero_weights = {}  # track_id -> position -> binding weight
    expectations = {}  # track_id -> expectation
    variances = {}  # track_id -> variance
    covariances = {}  # track_id -> other_id -> variance

    # Keep track of all intervals and weights for each track:
    orig_trackfile_handle = gzip.open(original_file) if original_file.endswith('gz') else open(original_file)
    header = None
    for track_line in orig_trackfile_handle:
        if track_line.startswith('#'):
            continue
        elif not header:
            header = track_line
            continue

        protein_id, track_id, _, current_intervals, current_nonzero_weights = track_line[:-1].split('\t')[:5]

        track_id = int(track_id)
        try:
            nonzero_weights[track_id] = {int(a.split(':')[0]): float(a.split(':')[1])
                                         for a in current_nonzero_weights.split(',')}
            intervals[track_id] = [(int(a.split('-')[0]), int(a.split('-')[1])) for a in current_intervals.split(',')]
        except IndexError:
            continue

    orig_trackfile_handle.close()

    # For each track, compute the expectation and variance:
    for track_id, track_intervals in sorted(intervals.items()):

        # current positions that we are interested in storing information for:
        current_positions = set()
        for interval_start, interval_end in track_intervals:
            for i in xrange(interval_start, interval_end + 1):
                current_positions.add(i)

        weight_vector = {i: nonzero_weights[track_id].get(i, 0) for i in current_positions}
        total_lambda = float(sum([mutation_rates.get(i, 1) for i in current_positions]))

        expectations[track_id] = sum(
            [mutation_rates.get(i, 1) * weight_vector[i] for i in current_positions]) / total_lambda
        variances[track_id] = sum([mutation_rates.get(i, 1) * (weight_vector[i] ** 2)
                                   for i in current_positions]) / total_lambda - expectations[track_id] ** 2

        # Compute COVARIANCE with each other track
        covariances[track_id] = {}
        for other_id, other_intervals in sorted(intervals.items()):
            if other_id <= track_id:
                continue

            # "other" positions
            other_positions = set()
            for interval_start, interval_end in other_intervals:
                for i in xrange(interval_start, interval_end + 1):
                    other_positions.add(i)

            # if no overlap, then covariance is 0, and move on
            overlap_positions = other_positions.intersection(current_positions)
            if len(overlap_positions) == 0:
                covariances[track_id][other_id] = '0.0|0.0|0.0'
                continue
            overlap_total_lambda = float(sum([mutation_rates.get(i, 1) for i in overlap_positions]))

            # Compute the covariance (over the overlapping section) as:
            #   (mutation likelihood) * SUM (w_i - E[w_i]) * (v_i - E[v_i])
            # where (1) E[w_i] and E[v_i] are calculated only over the overlapping section,

            other_weight_vector = {i: nonzero_weights[other_id].get(i, 0) for i in overlap_positions}
            other_total_lambda = float(sum([mutation_rates.get(i, 1) for i in other_positions]))

            likelihood_in_overlap = sum([mutation_rates.get(i, 1) for i in overlap_positions])

            likelihoods = (likelihood_in_overlap / total_lambda, likelihood_in_overlap / other_total_lambda)

            # if variance within the overlap region is 0, then covariance is 0, and move on
            if len(set([weight_vector[i] for i in overlap_positions])) < 2 or \
               len(set([other_weight_vector[i] for i in overlap_positions])) < 2:
                covariances[track_id][other_id] = '0.0' + '|' + str(likelihoods[0]) + '|' + str(likelihoods[1])
                continue

            mu_track_id = sum([mutation_rates.get(i, 1) * weight_vector[i] for i in
                               overlap_positions]) / (overlap_total_lambda if overlap_total_lambda > 0 else 1)
            mu_other_id = sum([mutation_rates.get(i, 1) * other_weight_vector[i] for i in
                               overlap_positions]) / (overlap_total_lambda if overlap_total_lambda > 0 else 1)

            overlap_covariance = sum([mutation_rates.get(i, 1) *
                                      (weight_vector[i] - mu_track_id) * (other_weight_vector[i] - mu_other_id)
                                      for i in overlap_positions]) / (
                                     overlap_total_lambda if overlap_total_lambda > 0 else 1)

            covariances[track_id][other_id] = '|'.join(map(str, [overlap_covariance, likelihoods[0], likelihoods[1]]))

    # Augment original file to include the covariances column
    out_handle = gzip.open(complete_file, 'w') if complete_file.endswith('gz') else open(complete_file, 'w')
    orig_trackfile_handle = gzip.open(original_file) if original_file.endswith('gz') else open(original_file)
    header = None
    for track_line in orig_trackfile_handle:
        if track_line.startswith('#'):
            out_handle.write(track_line)
            continue
        elif not header:
            header = track_line
            out_handle.write(header)
            continue

        track_id = int(track_line[:-1].split('\t')[1])
        try:
            out_handle.write('\t'.join(track_line[:-1].split('\t')[:5] +
                                       [str(expectations[track_id]),
                                        str(variances[track_id]),
                                        ','.join([str(other_id) + ':' + str(covariances[track_id][other_id])
                                                  for other_id in sorted(covariances[track_id].keys())])]) + '\n')
        except KeyError:
            continue
    orig_trackfile_handle.close()
    out_handle.close()


####################################################################################################
# STEP 3: Empirically determine the minimum number of mutations required to assume normality
#         so that we can use the Central Limit Theorem
####################################################################################################

def empirically_determine_minmuts(original_file, complete_file, mutation_rates=None):
    """
    :param original_file: full path to a tab-delimited file with (at least) columns [0] protein ID, [1] track ID,
                          [2] track name, [3] 0-index enrichment interval(s), [4] 0-index non-zero binding weights
    :param complete_file: full path to a new tab-delimited file with additional columns [5] expected value assuming
                          1 mutation, [6] variance and [7] covariances with all other tracks
    :param mutation_rates: list of protein length with values corresponding to likelihood of mutation
    :return: compute the expectation, variance, and covariances between all tracks and write out to new file
    """

    if not mutation_rates:  # index -> likelihood of mutation
        mutation_rates = {}  # no entries

    # Augment original file to include the minimum mutations column
    out_handle = gzip.open(complete_file, 'w') if complete_file.endswith('gz') else open(complete_file, 'w')
    orig_trackfile_handle = gzip.open(original_file) if original_file.endswith('gz') else open(original_file)
    header = None
    for track_line in orig_trackfile_handle:
        if track_line.startswith('#'):
            out_handle.write(track_line)
            continue
        elif not header:
            header = track_line
            out_handle.write(header)
            continue

        (protein_id, track_id, _, current_intervals, current_nonzero_weights,
         expectation, variance) = track_line[:-1].split('\t')[:7]

        try:
            nonzero_weights = {int(a.split(':')[0]): float(a.split(':')[1])
                               for a in current_nonzero_weights.split(',')}
            track_intervals = [(int(a.split('-')[0]), int(a.split('-')[1])) for a in current_intervals.split(',')]
            expectation = float(expectation)
            variance = float(variance)

        except IndexError:
            out_handle.write('\t'.join(track_line.strip().split('\t')[:7] +
                                       ['' if len(track_line.strip().split('\t')) < 8 else
                                        track_line.strip().split('\t')[7]]) + '\t500\t1\n')
            continue

        # Determine the minimum number of mutations for normally distributed z-scores

        # current positions that we are interested in storing information for:
        current_positions = set()
        for interval_start, interval_end in track_intervals:
            for i in xrange(interval_start, interval_end + 1):
                current_positions.add(i)

        weight_vector = {i: nonzero_weights.get(i, 0) for i in current_positions}
        total_lambda = float(sum([mutation_rates.get(i, 1) for i in current_positions]))

        if not variance > 0:
            variance = sum([mutation_rates.get(i, 1) * (weight_vector[i] ** 2)
                            for i in current_positions]) / total_lambda - expectation ** 2
        if not variance > 0:
            out_handle.write('\t'.join(track_line.strip().split('\t')[:7] +
                                       ['' if len(track_line.strip().split('\t')) < 8 else
                                        track_line.strip().split('\t')[7]]) + '\t500\t1\n')
            continue

        list_of_candidates = sorted(list(current_positions))  # all possible positions a mutation can land
        probability_distribution = [mutation_rates.get(i, 1) / total_lambda for i in list_of_candidates]

        minimum_muts = 500
        runtime_total = 1
        for mut_count in xrange(1, 501):  # test up to 500 mutations (but break as soon as we reach normality
            observed_zscores = []
            runtime_start = time.time()  # start the clock to measure performance
            for iteration in xrange(100):  # randomly sprinkle mutations 500 times then test for normality
                observed_zscores.append((sum([weight_vector[i] for i in np.random.choice(list_of_candidates, mut_count,
                                                                                         p=probability_distribution)]) -
                                         expectation) / variance)
            runtime_end = time.time() - runtime_start  # end the clock

            # now, are these z-scores normally distributed?
            if len(set(observed_zscores)) > 1 and shapiro(observed_zscores)[1] > 0.0005:
                minimum_muts = mut_count
                runtime_total = runtime_end
                break

        # and write it out...
        try:
            out_handle.write('\t'.join(track_line.strip().split('\t')[:7] +
                                       ['' if len(track_line.strip().split('\t')) < 8 else
                                        track_line.strip().split('\t')[7]]) +
                             '\t' + str(minimum_muts) + '\t' + str(runtime_total) + '\n')
        except KeyError:
            out_handle.write('\t'.join(track_line.strip().split('\t')[:7] +
                                       ['' if len(track_line.strip().split('\t')) < 8 else
                                        track_line.strip().split('\t')[7]]) + '\t500\t1\n')
            continue

    orig_trackfile_handle.close()
    out_handle.close()


####################################################################################################
# MAIN
####################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Create "weight vector" files for each protein.')

    parser.add_argument('--protein_list', type=str, default=data_path+'pertinint/modelable-proteins.txt',
                        help='Full path to a list of proteins (and their locations) on which to run.')
    parser.add_argument('--fasta_file', type=str, help='Full path to a fasta file with protein sequences (for length)',
                        default=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.pep.all.fa')
    parser.add_argument('--error_list', type=str, default=data_path+'pertinint/failed-modelable-proteins.txt',
                        help='Full path to a list of proteins (and their locations) which should be rerun.')
    parser.add_argument('--track_weights_path', type=str, default=data_path+'pertinint/track_weights/',
                        help='Full path to a directory where processed track weight files are stored, ' +
                             'separated into subdirectories by chromosome and gene.')
    parser.add_argument('--start', type=int, help='Starting protein index', default=0)
    parser.add_argument('--end', type=int, help='Ending protein index', default=None)
    parser.add_argument('--timeout', type=int, help='Number of seconds allowed for each protein to process',
                        default=600)

    parser.add_argument('--find_failed', dest='find_failed', action='store_true', default=False,
                        help='Get list of protein IDs that should have a track weight file but do not yet.')
    parser.add_argument('--create_manifest', dest='create_manifest', action='store_true', default=False,
                        help='Get list of all track weight files that have been successfully created.')
    parser.add_argument('--compress_files', dest='compress_files', action='store_true', default=False,
                        help='Compress track weight files using gzip (if not already done)')

    parser.add_argument('--initialize_only', dest='initialize_only', action='store_true', default=False)
    parser.add_argument('--significance_only', dest='significance_only', action='store_true', default=False)
    parser.add_argument('--minimum_muts_only', dest='minimum_muts_only', action='store_true', default=False)

    parser.add_argument('--empirical_test', dest='empirical_test', action='store_true', default=False)

    args = parser.parse_args()

    # ----------------------------------------------------------------------------------------------------
    # create PertInInt input/output directories as needed:
    if not os.path.isdir(data_path+'pertinint'):
        call(['mkdir', data_path+'pertinint'])
    if not args.track_weights_path.endswith('/'):
        args.track_weights_path += '/'
    if not os.path.isdir(args.track_weights_path):
        call(['mkdir', args.track_weights_path])

    # ----------------------------------------------------------------------------------------------------
    # confirm that FASTA file (needed in all steps) is readable
    if not os.path.isfile(args.fasta_file):
        sys.stderr.write('Could not find protein sequences in: '+args.fasta_file+'. Please run: \n' +
                         'python '+sys.argv[0]+' --fasta_file <full path to fasta-formatted sequence file>\n')
        sys.exit(1)

    # ----------------------------------------------------------------------------------------------------
    # (1) get the complete list of proteins to run on
    if not os.path.isfile(args.protein_list):
        sys.stderr.write('Could not find protein list: '+args.protein_list+'. Creating now... ')
        start = time.time()
        protein_list_file = find_modelable_prots([f[0] for f in standard_tracks] + binary_tracks + whole_gene_tracks,
                                                 args.fasta_file,
                                                 args.protein_list)
        sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n')

        if not protein_list_file:
            sys.stderr.write('Failed to create protein list!\n')
            sys.exit(1)

    sys.stderr.write('Reading in modelable proteins... ')
    start = time.time()
    proteins_to_process = read_modelable_prots(args.protein_list)
    sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n' +
                     '   > ' + args.protein_list + '\n\n')

    # ----------------------------------------------------------------------------------------------------
    # (2) find the list of failed protein IDs if specified to do so
    if args.find_failed:
        sys.stderr.write('Finding protein IDs that we should recreate track weight files for... ')
        start = time.time()
        if not os.path.isfile(args.protein_list):
            sys.stderr.write('Could not find protein list: ')
        error_file = find_failed_track_weights(args.protein_list,
                                               args.track_weights_path,
                                               args.error_list,
                                               None,  # set of chromosomes to restrict to
                                               long_proteins)  # ignore too-long genes without known cancer roles

        sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n' +
                         'Rerun failed track weight files by running: \n' +
                         'python '+sys.argv[0]+' --protein_list ' + error_file+'\n')
        sys.exit(0)

    # ----------------------------------------------------------------------------------------------------
    # (3) create a manifest file of existing track weight files if specified to do so
    if args.compress_files:
        sys.stderr.write('Gzipping track weight files where possible to save space... ')
        start = time.time()
        wd = os.getcwd()
        for chrom_dir in sorted([d for d in os.listdir(args.track_weights_path)
                                 if os.path.isdir(args.track_weights_path+d)]):
            for gene_dir in sorted([d for d in os.listdir(args.track_weights_path+chrom_dir)
                                    if os.path.isdir(args.track_weights_path+chrom_dir+'/'+d)]):
                for track_wt in sorted([d for d in os.listdir(args.track_weights_path+chrom_dir+'/'+gene_dir) if
                                        d.endswith('.trackweights.tsv')]):
                    os.chdir(args.track_weights_path+chrom_dir+'/'+gene_dir)
                    call(['gzip', track_wt])
        os.chdir(wd)
        sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n')
        sys.exit(0)

    # ----------------------------------------------------------------------------------------------------
    # (4) create a manifest file of existing track weight files if specified to do so
    if args.create_manifest:
        sys.stderr.write('Creating manifest file of all existing track weight files... ')
        start = time.time()

        # manifest files store all file names (one version considers only chromosomes 1-22, X, Y, MT
        manifest = open(args.track_weights_path+'manifest-all.log', 'w')
        manifest_trunc = open(args.track_weights_path+'manifest.log', 'w')

        # tarball files are used to create a tarball to move all data around
        pert_version = '0' if args.track_weights_path.endswith('track_weights/') else \
                       args.track_weights_path.split('/')[-2].replace('track_weights_', '')
        tarball_files = open('/'.join(args.track_weights_path.split('/')[:-2]) +
                             '/PertInInt-tracks_v'+pert_version+'-all.txt', 'w')
        tarball_files_trunc = open('/'.join(args.track_weights_path.split('/')[:-2]) +
                                   '/PertInInt-tracks_v'+pert_version+'.txt', 'w')

        dir_name = args.track_weights_path.split('/')[-2]+'/'
        tarball_files.write(dir_name + 'manifest.log\n')
        tarball_files.write(dir_name + 'manifest-all.log\n')
        tarball_files_trunc.write(dir_name + 'manifest.log\n')

        for chrom_dir in sorted([d for d in os.listdir(args.track_weights_path)
                                 if os.path.isdir(args.track_weights_path+d)]):
            for gene_dir in sorted([d for d in os.listdir(args.track_weights_path+chrom_dir)
                                    if os.path.isdir(args.track_weights_path+chrom_dir+'/'+d)]):
                for track_wt in sorted([d for d in os.listdir(args.track_weights_path+chrom_dir+'/'+gene_dir) if
                                        '.trackweights.tsv' in d]):
                    manifest.write(chrom_dir+'/'+gene_dir+'/'+track_wt+'\n')
                    tarball_files.write(dir_name+chrom_dir+'/'+gene_dir+'/'+track_wt+'\n')
                    if chrom_dir in map(str, range(1, 23)) + ['X', 'Y']:
                        manifest_trunc.write(chrom_dir+'/'+gene_dir+'/'+track_wt+'\n')
                        tarball_files_trunc.write(dir_name+chrom_dir+'/'+gene_dir+'/'+track_wt+'\n')
        manifest.close()
        tarball_files.close()
        tarball_files_trunc.close()

        sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n' +
                         '   > '+args.track_weights_path + 'manifest.log\n' +
                         '   > ' + args.track_weights_path + 'manifest-all.log\n' +
                         '   > create_tarball$ cd '+'/'.join(args.track_weights_path.split('/')[:-2])+'/\n' +
                         '   > create_tarball$ tar -cvzf PertInInt-tracks_v0.tar.gz -T PertInInt-tracks_v0.txt\n' +
                         '   > create_tarball$ tar -cvzf PertInInt-tracks_v0-all.tar.gz -T PertInInt-tracks_v0-all.txt\n\n')
        sys.exit(0)

    # ----------------------------------------------------------------------------------------------------
    # (5) restrict proteins to run on if specified
    start_index = min(max(0, args.start), len(proteins_to_process) - 1) if args.start else 0
    end_index = min(args.end, len(proteins_to_process)) if args.end else len(proteins_to_process)
    proteins_to_process = proteins_to_process[start_index:end_index]

    # drop TTN, OBSCN, NEB because they take WAY too long
    proteins_to_process = [prot for prot in proteins_to_process if not prot.split('/')[1] in long_proteins]

    # ----------------------------------------------------------------------------------------------------
    # TEST CASE:
    # proteins_to_process = ['12/ENSG00000135390/ENSP00000377878']
    # proteins_to_process = ['13/ENSG00000102554/ENSP00000366915']  # KLF5
    # ----------------------------------------------------------------------------------------------------

    if len(proteins_to_process) < 1:
        sys.exit(0)

    # ----------------------------------------------------------------------------------------------------
    # (6) initialize all track weight files (i.e., names, intervals, and functional weights)
    if not (args.significance_only or args.minimum_muts_only):

        # (5.1) read in sequence lengths for all proteins to be processed
        sys.stderr.write('Loading sequence lengths... ')
        start = time.time()
        sequence_lengths = sequence_lens(args.fasta_file, [prot.split('/')[2] for prot in proteins_to_process])
        sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n' +
                         '   > ' + args.fasta_file + '\n\n')

        # (5.2) read in all standard, binary, and whole gene track information
        sys.stderr.write('Loading standard, binary, and whole gene tracks... ')
        start = time.time()
        standard_regions, region_weights = process_standard_tracks(standard_tracks)
        binary_regions = process_binary_tracks(binary_tracks)
        gene_count, gene_weights = process_whole_gene_tracks(whole_gene_tracks)
        sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n' +
                         '\n'.join(['   > ' + filename for filename in
                                    [item for sublist in standard_tracks for item in sublist] +
                                    binary_tracks + whole_gene_tracks]) + '\n\n')

        # (5.3) for each protein, create a new "track weights" file:
        for seq_index, full_prot_loc in enumerate(proteins_to_process):
            chromosome, gene_id, prot_id = full_prot_loc.split('/')

            signal.signal(signal.SIGALRM, handler)  # Register the signal function handler
            signal.alarm(args.timeout)  # Define a timeout for this function

            step1_file = args.track_weights_path + full_prot_loc + '.trackweights-orig.tsv'

            try:
                sys.stderr.write('['+str(seq_index)+'] Creating track weight file for ' +
                                 full_prot_loc + '.trackweights-orig.tsv... ')
                start = time.time()
                process_domain_weights(prot_id,  # protein ID
                                       gene_id,  # gene ID
                                       chromosome,  # chromosome
                                       sequence_lengths[prot_id],  # protein length
                                       step1_file,  # full path to output file
                                       standard_regions.get(prot_id, None),  # standard track locations
                                       region_weights,  # standard track functional scores
                                       binary_regions.get(prot_id, None),  # binary track locations
                                       gene_count,  # total genes evaluated
                                       gene_weights.get(prot_id, None))  # whole gene relative mutability
                signal.alarm(0)  # Cancel the alarm if we made it to this point
                sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n')

            except Exception as exc:
                sys.stderr.write(str(exc) + ' on ' + full_prot_loc + '\n')
                continue

    # ----------------------------------------------------------------------------------------------------
    # (7) Model the likelihood of per-position mutations according to some background model
    all_mutation_rates = None
    if not args.initialize_only:
        sys.stderr.write('Loading protein position -> mutation likelihood mappings... ')
        start = time.time()
        all_mutation_rates = mutational_likelihoods_overall(mutation_biases, proteins_to_process)
        sys.stderr.write('finished in '+reformat_time(time.time()-start)+'!\n' +
                         '   > ' + mutation_biases + '\n\n')

    # ----------------------------------------------------------------------------------------------------
    # (8) include expectation, variance, and covariance estimates in the track weight files
    if not (args.initialize_only or args.minimum_muts_only):

        for seq_index, full_prot_loc in enumerate(proteins_to_process):
            chromosome, gene_id, prot_id = full_prot_loc.split('/')

            signal.signal(signal.SIGALRM, handler)  # Register the signal function handler
            signal.alarm(args.timeout)  # Define a timeout for this function

            step1_file = args.track_weights_path + full_prot_loc + '.trackweights-orig.tsv'

            if not os.path.isfile(step1_file):
                sys.stderr.write('Could not find '+step1_file+'\n')
                continue

            step2_file = args.track_weights_path + full_prot_loc + '.trackweights-significance.tsv'

            try:
                sys.stderr.write('['+str(seq_index)+'] Computing significance in track file ' +
                                 full_prot_loc + '.trackweights-significance.tsv... ')
                start = time.time()
                # (7.1) convert background mutational likelihoods into per-protein (binned) values:
                current_mutation_rates = {index: variation for index, variation in
                                          enumerate(mutational_likelihoods_by_protein(all_mutation_rates,
                                                                                      prot_id,
                                                                                      args.fasta_file))}

                # (7.2) include expectation, variance, and covariance in original files:
                update_domain_weights(step1_file, step2_file, current_mutation_rates)

                # (7.3) save space by removing original file (if we successfully reached this point):
                if os.path.isfile(step1_file) and os.path.isfile(step2_file) and step1_file != step2_file:
                    call(['rm', step1_file])

                signal.alarm(0)  # Cancel the alarm if we made it to this point
                sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n')

            except Exception as exc:
                sys.stderr.write(str(exc) + ' on ' + full_prot_loc + '\n')
                continue

    # ----------------------------------------------------------------------------------------------------
    # (9) empirically determine the minimum number of mutations required for each track to assume normality
    if not (args.initialize_only or args.significance_only):

        for seq_index, full_prot_loc in enumerate(proteins_to_process):
            chromosome, gene_id, prot_id = full_prot_loc.split('/')

            signal.signal(signal.SIGALRM, handler)  # Register the signal function handler
            signal.alarm(args.timeout)  # Define a timeout for this function

            step2_file = args.track_weights_path + full_prot_loc + '.trackweights-significance.tsv'

            if not os.path.isfile(step2_file):
                sys.stderr.write('Could not find '+step2_file+'\n')
                continue

            step3_file = args.track_weights_path + full_prot_loc + '.trackweights.tsv'

            try:
                sys.stderr.write('['+str(seq_index)+'] Computing normality for track file ' +
                                 full_prot_loc + '.trackweights.tsv... ')
                start = time.time()
                # (8.1) convert background mutational likelihoods into per-protein (binned) values:
                current_mutation_rates = {index: variation for index, variation in
                                          enumerate(mutational_likelihoods_by_protein(all_mutation_rates,
                                                                                      prot_id,
                                                                                      args.fasta_file))}

                # (8.2) empirically determine the minimum number of mutations required to assume normality:
                empirically_determine_minmuts(step2_file, step3_file, current_mutation_rates)

                # (8.3) save space by removing original file (if we successfully reached this point):
                if os.path.isfile(step2_file) and os.path.isfile(step3_file) and step2_file != step3_file:
                    call(['rm', step2_file])

                signal.alarm(0)  # Cancel the alarm if we made it to this point
                sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n')

            except Exception as exc:
                sys.stderr.write(str(exc) + ' on ' + full_prot_loc + '\n')
                continue
