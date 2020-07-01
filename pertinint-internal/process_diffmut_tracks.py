#!/usr/bin/python

"""
Download and process whole-gene tracks based on 1000 Genome variation binned and processed per-person
by DiffMut:
Przytycki, P.F. and Singh, M. (2017). "Differential analysis between somatic mutation and germline variation
profiles reveals cancer-related genes." Genome Med, 9: 79. doi: 10.1186/s13073-017-0465-6

Contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import math
import time
import gzip
import argparse
from config import data_path
from subprocess import call


####################################################################################################
# DOWNLOAD AND PROCESS DIFFMUT INPUT DATA
####################################################################################################

def reformat_time(run_time):
    """
    :param run_time: total time elapsed in seconds
    :return: string corresponding to a properly formatted (days, hours, minutes, seconds) time
    """

    m, s = divmod(run_time, 60)
    h, m = divmod(m, 60)
    d, h = divmod(h, 24)
    return ':'.join(map(lambda x: str(int(x)).zfill(2), [d, h, m, s]))


####################################################################################################

def download_diffmut_rawdata(outdir):
    """
    :param outdir: full path to an output directory to store two downloaded files (from DiffMut)
    :return: full paths to the two downloaded files
    """

    data_matrix = outdir + 'diffmut_binned_natural_variation.txt'
    data_labels = outdir + 'diffmut_protein_names.txt'

    # R-readable matrix of 19,460 genes x 100 processed bins
    if not os.path.isfile(data_matrix):
        call(['wget',
              'https://raw.githubusercontent.com/Singh-Lab/Differential-Mutation-Analysis/master/Data/' +
              'natDistBinned.txt',
              '-O',
              data_matrix])

    # corresponding list of ordered protein names
    if not os.path.isfile(data_labels):
        call(['wget',
              'https://raw.githubusercontent.com/Singh-Lab/Differential-Mutation-Analysis/master/Data/protNames.txt',
              '-O',
              data_labels])

    if not os.path.isfile(data_matrix):
        sys.stderr.write('Failed to download data matrix from https://github.com/Singh-Lab/' +
                         'Differential-Mutation-Analysis... exiting\n')
        sys.exit(1)

    if not os.path.isfile(data_labels):
        sys.stderr.write('Failed to download data labels from https://github.com/Singh-Lab/' +
                         'Differential-Mutation-Analysis... exiting\n')
        sys.exit(1)

    return data_matrix, data_labels


####################################################################################################

def protein_list(protnames_file):
    """
    :param protnames_file: full path to the protein file ordered for use by DiffMut
    :return: the list of proteins specified in the file provided by DiffMut (parsed from 1000 Genomes)
    """

    prot_handle = gzip.open(protnames_file) if protnames_file.endswith('gz') else open(protnames_file)
    prot_names = [prot_line.strip() for prot_line in prot_handle]
    prot_handle.close()

    return prot_names


####################################################################################################

def gene_probabilities_vector(binnedranks, protlist):
    """
    :param binnedranks: full path to a file containing the binned ranks of genes
    :return: dictionary with keys corresponding to gene names and values corresponding to that gene's
             likelihood of harboring a missense mutation
    """

    # get the ordered list of genes:
    ordered_genes = protein_list(protlist)
    gene_to_burden = {gene: 0 for gene in ordered_genes}  # keep track of the importance of each gene

    # get the relative importance of the gene. First row is least important, last row is more important
    burden_handle = gzip.open(binnedranks) if binnedranks.endswith('gz') else open(binnedranks)
    current_row = 1
    for row in burden_handle:
        burdens = row[:-1].split()  # this length should match the number of genes
        for i in xrange(len(burdens)):
            gene_to_burden[ordered_genes[i]] += float(burdens[i]) * current_row
        current_row += 1  # last row means top 1% of mutations
    burden_handle.close()

    # normalize the dictionary for per-gene probabilities
    total_burden = sum(gene_to_burden.values())
    for gene_id in gene_to_burden.keys():
        gene_to_burden[gene_id] = gene_to_burden[gene_id] / total_burden

    return gene_to_burden


####################################################################################################
# CREATE WHOLE GENE TRACK OUTPUT FILE
####################################################################################################

def get_gene_aliases(gene_info_file=data_path+'diffmut/ncbi_gene_info.txt.gz'):
    """
    :param gene_info_file: full path to the diffmut data directory where a gene alias file will be
                        downloaded (if needed) and parsed
    :return: dictionary of hugo gene name -> set(all other aliases)
    """

    if not os.path.isfile(gene_info_file):
        call(['wget',
              'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz',
              '-O',
              gene_info_file])

    hugo_mapping = {}

    gene_info_handle = gzip.open(gene_info_file, 'rt') if gene_info_file.endswith('gz') else open(gene_info_file)
    for gene_line in gene_info_handle:
        if gene_line.startswith('#'):
            continue

        species, _, symbol, _, synonyms = gene_line[:-1].split('\t')[:5]
        if species == '9606':
            all_names = set([symbol] + synonyms.split('|')) if synonyms != '-' else {symbol}

            for hugo_id in all_names:
                if hugo_id not in hugo_mapping:
                    hugo_mapping[hugo_id] = set()
                hugo_mapping[hugo_id].add(symbol)  # = hugo_mapping[hugo_id].union(all_names)

    gene_info_handle.close()

    return hugo_mapping


####################################################################################################

def get_protein_mapping(humanprots):
    """
    :param humanprots: full path to a tab-delimited file where Ensembl gene IDs, protein IDs and HGNC symbols are
    :return: two dictionaries of hugo ID -> set(ensembl gene IDs) and ensembl gene ID -> set(ensembl prot IDs)
    """

    hugo_mapping = {}
    gene_mapping = {}

    mapping_handle = gzip.open(humanprots) if humanprots.endswith('gz') else open(humanprots)
    header = None
    for fline in mapping_handle:
        if fline.startswith('#'):
            continue
        elif not header:
            header = fline[:-1].split('\t')
            continue

        v = fline[:-1].split('\t')
        prot_id = v[header.index('Protein stable ID')]
        gene_id = v[header.index('Gene stable ID')]
        hugo_ids = set([nm for nm in [v[header.index('HGNC symbol')],
                                      v[header.index('Gene name')]] if len(nm.strip()) > 0])

        for hugo_id in hugo_ids:
            if hugo_id not in hugo_mapping:
                hugo_mapping[hugo_id] = set()
            hugo_mapping[hugo_id].add(gene_id)

            if gene_id not in gene_mapping:
                gene_mapping[gene_id] = set()
            gene_mapping[gene_id].add(prot_id)
    mapping_handle.close()

    return hugo_mapping, gene_mapping


####################################################################################################

def output_per_prot_likelihoods(binnedranks_file, protnames_file, mapping_file, gene_info_file, out_track_file,
                                log_file):
    """
    :param binnedranks_file: full path to a file containing the binned ranks of genes
    :param protnames_file: full path to the protein file ordered for use by DiffMut
    :param mapping_file: full path to a tab-delimited file with Ensembl gene IDs, protein IDs, and HGNC IDs
    :param gene_info_file: full path to a file containing all gene name synonyms
    :param out_track_file: full path to file to write to
    :return: none but write appropriate information to file
    """

    # columns will be Ensembl protein ID, Ensembl gene ID, total genes, likelihood, description (hugo)
    # (1) get the mappings from Hugo ID -> Ensembl Gene IDs and from Ensembl Gene IDs -> Ensembl Protein IDs:
    # (3) get Hugo ID synonyms for proper mapping
    sys.stderr.write('Obtaining HGNC symbol aliases for all human genes... ')
    time_start = time.time()
    hugo_mapping, gene_mapping = get_protein_mapping(mapping_file)
    hugo_aliases = get_gene_aliases(gene_info_file)
    sys.stderr.write('finished in '+reformat_time(time.time()-time_start)+'!\n')

    # (2) get the per-gene relative likelihoods of mutation:
    hugo_to_mutlikelihood = gene_probabilities_vector(binnedranks_file, protnames_file)

    # (3) open the output file to write to
    out_handle = gzip.open(out_track_file, 'wt') if out_track_file.endswith('gz') else open(out_track_file, 'w')
    out_handle.write('\n'.join(['# Per-gene likelihood of harboring a missense mutation derived from 1000 Genomes data',
                                '# Preprocessed by DiffMut (Przytycki & Singh 2017). Original files include:',
                                '# ' + binnedranks_file,
                                '# ' + protnames_file,
                                '# Gene name mappings found in:',
                                '# ' + mapping_file,
                                '# ' + gene_info_file,
                                '\t'.join(['#Ensembl_Protein_ID', 'Ensembl_Gene_ID', 'Total_HGNC_Genes_Evaluated',
                                           'Relative_Mutational_Likelihood', 'Description'])])+'\n')

    log_handle = gzip.open(log_file, 'wt') if log_file.endswith('gz') else open(log_file, 'w')

    # (4) get total number of genes
    ensembl_gene_set = set()  # total Ensembl genes with associated values
    hgnc_attempted = set()  # total HGNC gene symbols evaluated by DiffMut
    hgnc_found = set()  # total HGNC gene symbols that are still protein-coding genes (and could be mapped)
    for hugo_id in hugo_to_mutlikelihood.keys():
        hgnc_attempted.add(hugo_id)

        # if there is a clean, 1-to-1 mapping to Ensembl gene ID:
        if hugo_id in hugo_mapping:
            hgnc_found.add(hugo_id)
            ensembl_gene_set = ensembl_gene_set.union(hugo_mapping[hugo_id])

        # otherwise, if we could find the "primary" gene name for this hugo ID:
        elif hugo_id in hugo_aliases:
            for hugo_alias in hugo_aliases[hugo_id]:
                if hugo_alias in hugo_mapping and hugo_alias not in hugo_to_mutlikelihood.keys():
                    hgnc_found.add(hugo_alias)
                    ensembl_gene_set = ensembl_gene_set.union(hugo_mapping[hugo_alias])

    # (5) write out information
    for hugo_id, mut_likelihood in sorted(hugo_to_mutlikelihood.items()):

        # get all the aliases for this hugo_id:
        if hugo_id in hugo_mapping:
            final_id = hugo_id
            gene_set = hugo_mapping[hugo_id]
        else:
            gene_set = set()
            final_ids = set()
            for hugo_alias in hugo_aliases.get(hugo_id, set()):
                if hugo_alias in hugo_mapping and hugo_alias not in hugo_to_mutlikelihood.keys():
                    final_ids.add(hugo_alias)
                    gene_set = gene_set.union(hugo_mapping[hugo_alias])
            final_id = '|'.join(sorted(list(final_ids)))

        gene_set = sorted(list(gene_set))
        if len(gene_set) < 1:
            log_handle.write('|'.join(sorted(list(hugo_aliases.get(hugo_id, {hugo_id}))))+'\n')

        for gene_id in gene_set:
            for prot_id in sorted(list(gene_mapping.get(gene_id, set()))):
                if len(prot_id.strip()) > 0:
                    out_handle.write('\t'.join([prot_id,
                                                gene_id,
                                                str(len(hgnc_attempted)),
                                                str(mut_likelihood),
                                                'hgnc_symbol=' + final_id + ';' +
                                                'original_diffmut_symbol=' + hugo_id + ';' +
                                                'ensembl_genes_included=' + str(len(ensembl_gene_set)) + ';' +
                                                'hgnc_genes_included=' + str(len(hgnc_found))]) + '\n')
    out_handle.close()
    log_handle.close()

    return out_track_file, log_file


####################################################################################################
# PROCESS MUTATION DATA
####################################################################################################

def get_full_protids(humanprots):
    """
    :param humanprots: full path to a fasta-formatted file where ensembl protein IDs and Hugo gene names
                       are in the sequence headers
    :return: dictionary of truncated Ensembl protein IDs to full IDs (with versioning)
    """

    full_prot_ids = {}
    fasta_handle = gzip.open(humanprots) if humanprots.endswith('gz') else open(humanprots)
    for fasta_line in fasta_handle:
        if fasta_line.startswith('>'):
            prot_id = fasta_line[1:-1].split()[0]
            full_prot_ids[prot_id.split('.')[0]] = prot_id
    fasta_handle.close()

    return full_prot_ids


####################################################################################################

def mutcounts(mutfile, limit_expression=True, silent_flag=False, include_samples=None, expression_file=None,
              seq_file=None):
    """
    :param mutfile: full path to a mutation file (.maf format)
    :return: dictionary of gene -> [missense mutation count per individual]
    """

    expression_by_gene = {}
    if limit_expression:
        if not os.path.isfile(expression_file):
            sys.stderr.write('Attempted to limit to mutations in genes expressed at TPM > 0.1, ' +
                             'but could not find formatted expression file in ' + str(expression_file) + '\n')
        else:
            with gzip.open(expression_file) as exp_handle:
                for exp_line in exp_handle:
                    if exp_line.startswith('#'):
                        continue
                    gene_name = exp_line.split('\t')[0].split(',')[1]
                    expression_by_gene[gene_name] = set(exp_line.split('\t')[1].split(','))

    gene_names = protein_list()  # genes that we have natural variation information for
    tumor_samples = set()  # keep track of all tumor samples (for eventual sorted list with mut counts)

    mutation_counts = {gn: {} for gn in gene_names}  # gene -> sample_id -> mutation count
    mutation_values = {}  # prot -> (total mutation value)
    prots = get_full_protids(seq_file)
    total_mutations, total_mutational_value = 0, 0.

    # ------------------------------------------------------------------------------------------------
    header = None  # we need to know the gene name, mutation type, and tumor sample ID
    mut_handle = gzip.open(mutfile) if mutfile.endswith('gz') else open(mutfile)
    for mutline in mut_handle:
        if mutline.startswith('#'):
            continue
        if not header:
            header = mutline.lower()[:-1].split('\t')
            continue

        protid = mutline[:-1].split('\t')[header.index('ensp')]
        if protid not in prots:
            continue
        protid = prots[protid]  # reset to the versioning that we use

        gene_name = mutline[:-1].split('\t')[header.index('hugo_symbol')]
        mut_type = mutline[:-1].split('\t')[header.index('variant_classification')]
        sample_id = '-'.join(mutline[:-1].split('\t')[header.index('tumor_sample_barcode')].split('-')[:4])
        mut_val = float(mutline[:-1].split('\t')[header.index('t_alt_count')]) / \
                  float(mutline[:-1].split('\t')[header.index('t_depth')])

        try:
            aachange = mutline[:-1].split('\t')[header.index('hgvsp_short')][2:]  # e.g., p.L414L
            int(''.join([i for i in list(aachange) if i in map(str, range(10))])) - 1
        except ValueError:
            continue

        if include_samples and sample_id not in include_samples:
            continue

        if len(expression_by_gene.keys()) > 0 and sample_id not in expression_by_gene.get(gene_name, set()):
            continue

        if (silent_flag and 'silent' in mut_type.lower()) or \
            (not silent_flag and ('missense' in mut_type.lower() or 'nonsense' in mut_type.lower())):
            total_mutations += 1
            total_mutational_value += mut_val

            if gene_name in gene_names:
                if protid not in mutation_values:
                    mutation_values[protid] = [gene_name, 0.]
                mutation_values[protid][1] += mut_val

                tumor_samples.add(sample_id)
                if sample_id not in mutation_counts[gene_name]:
                    mutation_counts[gene_name][sample_id] = 0
                mutation_counts[gene_name][sample_id] += 1
    mut_handle.close()

    sorted_samples = sorted(list(tumor_samples))
    final_mut_counts = {gn: [mutation_counts[gn].get(sid, 0) for sid in sorted_samples] for gn in gene_names}

    return final_mut_counts, mutation_values, float(total_mutational_value), float(total_mutations)


####################################################################################################
# CALCULATE Z-SCORES
####################################################################################################

def calculate_original_zscores(gene_probabilities, gene_to_mutval, total_mut_count, total_mut_val):
    """
    :param gene_probabilities: dictionary of gene -> probability of mutation
    :param gene_to_mutval: dictionary of gene -> total mutational value
    :param total_mutations: total number of nonsynonymous mutations observed in file
    :return: dictionary of gene -> zscore
    """

    # determine how to scale the total mutation count to get reasonably-sized Z-scores

    exp_mutval = total_mut_val / total_mut_count  # average mutation value
    scaled_mutcount = math.sqrt(total_mut_val) / exp_mutval  # scale the total mutations down to its square root
    total_mutval = scaled_mutcount * exp_mutval  # reset the total mutational value
    totsq_mutval = scaled_mutcount * exp_mutval * exp_mutval  # reset the total mutational squared values
    scale_factor = math.sqrt(total_mut_val) / total_mut_val  # scale down the eventual observed mut counts

    gene_to_zscore = {}

    for prot_id, (gene_name, mut_observed) in gene_to_mutval.items():
        # whole gene track looks like: [0,0,0,0,...,1,...,0,0,0,0]. The sum of this AND the sum of these entries squared
        # are both 1, which is why the variance calculation works out below (the likelihood of landing on the 1 is same)
        mut_expected = gene_probabilities[gene_name]
        mut_variance = mut_expected - mut_expected ** 2

        prot_observed = mut_observed * scale_factor
        prot_expected = total_mutval * mut_expected
        prot_variance = totsq_mutval * mut_variance

        gene_to_zscore[prot_id] = (prot_observed - prot_expected) / math.sqrt(
            prot_variance if prot_variance > 0 else 1.)

    return gene_to_zscore


####################################################################################################

def whole_gene_zscores(maf_file,
                       diffmut_path=data_path+'diffmut/',
                       limit_expression=True,
                       silent_flag=False,
                       include_samples=None,
                       expression_file=data_path+'gdc/expression/TCGA_GRCh38_expressed-genes_TPM.tsv.gz',
                       seq_file=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.pep.all.withgenelocs.fa.gz'):
    """
    :param maf_file: string corresponding to the type of cancer we are expecting
    :param diffmut_path: full path to where reprocessed DiffMut input files are stored
    :param limit_expression: boolean indicating whether to restrict to "expressed" genes (True) or not
    :param silent_flag: boolean indicating whether to restrict to silent mutations (True) or not
    :param include_samples: set of samples to include (if None, then include all samples)
    :param expression_file: full path to a file containing information about expressed genes
    :param seq_file: full path to the FASTA file containing all genes we are modeling
    :return: dictionary of hugo gene symbol -> whole gene zscore
    """

    # try to process the mutation data:
    (gene_to_samplemuts,
     gene_to_mutval,
     total_mutval,
     total_mutations) = mutcounts(
        maf_file,
        limit_expression,
        silent_flag,
        include_samples,
        expression_file,
        seq_file
    )

    # get the probability of each gene being mutated:
    gene_probabilities = gene_probabilities_vector(diffmut_path+'diffmut_binned_natural_variation.txt',
                                                   diffmut_path+'diffmut_protein_names.txt')

    return calculate_original_zscores(gene_probabilities, gene_to_mutval, total_mutations, total_mutval)


####################################################################################################
# MAIN
####################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process DiffMut raw data to compute whole-gene Z-scores.')
    parser.add_argument('--diffmut_path', type=str, default=data_path + 'diffmut/',
                        help='Full path to a directory containing two raw files downloaded from DiffMut')
    parser.add_argument('--gene_mapping', type=str,
                        default=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.toHGNC.tsv',
                        help='Full path to a tab-delimited file containing a mapping from Ensembl gene and protein ' +
                             'IDs to HGNC gene symbols')
    parser.add_argument('--expression_file', type=str,
                        help='Full path to a tab-delimited file containing lists of genes that are expressed in ' +
                             'particular tumor samples',
                        default=data_path+'gdc/expression/TCGA_GRCh38_expressed-genes_TPM.tsv.gz')
    args = parser.parse_args()

    # (1) download data if it is not already found:
    sys.stderr.write('Downloading input data files from DiffMut GitHub... ')
    start = time.time()
    if not os.path.isdir(args.diffmut_path):
        call(['mkdir', args.diffmut_path])
    binned_ranks_file, gene_names_file = download_diffmut_rawdata(args.diffmut_path)
    sys.stderr.write('finished in '+reformat_time(time.time()-start)+'!\n')

    # (2) confirm that mapping file exists
    if not os.path.isfile(args.gene_mapping):
        sys.stderr.write('Could not find tab-delimited mapping file:\n'+args.gene_mapping+'...Exiting.\n')
        sys.exit(1)

    # (3) create output track file
    sys.stderr.write('Creating DiffMut natural variation track file...\n')
    start = time.time()
    (results_file,
     log_file) = output_per_prot_likelihoods(binned_ranks_file,
                                             gene_names_file,
                                             args.gene_mapping,
                                             args.diffmut_path+'ncbi_gene_info.txt.gz',
                                             args.diffmut_path+'diffmut-naturalvariation_wholegene-GRCh38.txt.gz',
                                             args.diffmut_path+'diffmut-naturalvariation_wholegene-GRCh38.' +
                                             'mapping-errors')
    sys.stderr.write('finished in ' + reformat_time(time.time() - start) + '!\n')

    # (5) print out locations of output files to inspect
    sys.stderr.write('Results in: '+results_file+'\n')
    sys.stderr.write('Errors in: '+log_file+'\n')
