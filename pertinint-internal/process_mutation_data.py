#!/usr/bin/python

"""
Obtain and parse variant .maf files to extract relevant mutations

Contact Shilpa Nadimpalli Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import gzip
import math
import argparse
from config import data_path, GENOME_BUILD as BUILD
from subprocess import call
from intervaltree import IntervalTree

########################################################################################################
# CONSTANTS
########################################################################################################

GENCODE = {
    'AGC': 'S', 'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TCN': 'S',  # Serine
    'AGA': 'R', 'AGG': 'R', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'CGN': 'R',  # Arginine
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CTN': 'L', 'TTA': 'L', 'TTG': 'L',  # Leucine
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'ACN': 'T',  # Threonine
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CCN': 'P',  # Proline
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GTN': 'V',  # Valine
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GCN': 'A',  # Alanine
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'GGN': 'G',  # Glycine
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I',  # Isoleucine
    'TAA': '_', 'TAG': '_', 'TGA': '_',  # Stop
    'AAC': 'N', 'AAT': 'N',  # Asparagine
    'AAA': 'K', 'AAG': 'K',  # Lysine
    'CAC': 'H', 'CAT': 'H',  # Histidine
    'CAA': 'Q', 'CAG': 'Q',  # Glutamine
    'GAC': 'D', 'GAT': 'D',  # Aspartic Acid
    'GAA': 'E', 'GAG': 'E',  # Glutamic Acid
    'TTC': 'F', 'TTT': 'F',  # Phenylalanine
    'TAC': 'Y', 'TAT': 'Y',  # Tyrosine
    'TGC': 'C', 'TGT': 'C',  # Cysteine
    'ATG': 'M',  # Methionine
    'TGG': 'W',  # Tryptophan
}

COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
              'a': 't', 't': 'a', 'c': 'g', 'g': 'c', 'n': 'n'}


########################################################################################################
# HELPER (MATH) FUNCTIONS
########################################################################################################

def is_numeric(s):
    """
    :param s: an arbitrary string
    :return: True if the string can be converted to a float and *is not NAN*, False otherwise
    """
    try:
        float_val = float(s)
        if str(float_val).lower() != 'nan':
            return True
        else:
            return False
    except ValueError:
        return False


########################################################################################################

def arithmetic_mean(value_list):
    """
    :param value_list: list of values
    :return: the arithmetic average of the list
    """

    all_numeric_values = [float(entry) for entry in value_list if is_numeric(entry)]
    return sum(all_numeric_values) / max(len(all_numeric_values), 1)


########################################################################################################

def standard_deviation(value_list):
    """
    :param value_list: list of values
    :return: the standard deviation of the values in the list
    """

    list_average = arithmetic_mean(value_list)
    return math.sqrt(
        sum((float(entry) - list_average) ** 2 for entry in value_list if is_numeric(entry)) / list_average)


########################################################################################################
# AGGREGATE DATA
########################################################################################################

def aggregate_groups():
    """
    :return: dictionary of aggregate cancer name -> set of cancer subtypes that are included
    """

    return {'Aggregate': {'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC',
                          'KICH', 'KIRC', 'KIRP', 'LAML', 'LCML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV',
                          'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM',
                          'UCEC', 'UCS', 'UVM'},
            'AggregateTrunc': {'ACC', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'GBM', 'HNSC',
                               'KICH', 'KIRC', 'KIRP', 'LAML', 'LCML', 'LGG', 'LIHC', 'MESO', 'OV',
                               'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'TGCT', 'THCA', 'THYM',
                               'UCEC', 'UCS', 'UVM'},
            'PANGI': {'ESCA', 'STAD', 'COAD', 'READ'},
            'KIPAN': {'KICH', 'KIRC', 'KIRP'},
            'COADREAD': {'COAD', 'READ'},
            'GBMLGG': {'GBM', 'LGG'},
            'STES': {'STAD', 'ESCA'}}


########################################################################################################

def create_aggregate_maf_files(input_mafs, output_maf):
    """
    :param input_mafs: set or list of full paths to maf files to be included in the output maf file
    :param output_maf: full path to an "output" maf file to write aggregate mutations to
    :return: none; this step might require a LOT of memory!
    """

    # (1) process JUST the headers from each of the input mafs to make sure we combine columns correctly
    comments = ['# All mutation data from the following files combined in ' + output_maf + ':']
    master_header = []
    for maf_file in input_mafs:
        if not os.path.isfile(maf_file):
            sys.stderr.write('No maf file ' + maf_file + '\n')
            continue

        current_header = None
        maf_handle = gzip.open(maf_file, 'rt') if maf_file.endswith('gz') else open(maf_file)
        for maf_line in maf_handle:
            if maf_line.startswith('#'):
                continue
            current_header = maf_line[:-1].split('\t')
            break
        maf_handle.close()

        # add all elements of the new header into the master_header
        new_elements = [elem for elem in current_header if elem not in master_header]
        master_header.extend(new_elements)

        # keep track of a running comment header for the output file
        comments.append('# ' + maf_file)

    # (2) NOW we are ready to combine all the maf files
    if len(master_header) < 1:
        sys.stderr.write('No input maf files to be written to ' + output_maf + '!\n')
        return None

    master_handle = gzip.open(output_maf, 'wt') if output_maf.endswith('gz') else open(output_maf, 'w')
    master_handle.write('\n'.join(comments) + '\n')
    master_handle.write('\t'.join(master_header) + '\n')

    seen = set()  # keep track of mutations that we have already seen (to avoid duplicates as we go)
    for maf_file in input_mafs:
        if not os.path.isfile(maf_file):
            continue

        current_header = None
        maf_handle = gzip.open(maf_file, 'rt') if maf_file.endswith('gz') else open(maf_file)
        for maf_line in maf_handle:
            if maf_line.startswith('#'):
                continue
            elif not current_header:
                current_header = maf_line[:-1].split('\t')
                continue
            cvals = maf_line[:-1].split('\t')
            ovals = [cvals[current_header.index(elem)] if elem in current_header else '' for elem in master_header]
            curr_mut = [cvals[current_header.index(elem)] for elem in
                        ['NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Tumor_Sample_Barcode',
                         'Reference_Allele']]
            curr_mut.append(','.join(sorted(list({cvals[current_header.index('Tumor_Seq_Allele1')],
                                                  cvals[current_header.index('Tumor_Seq_Allele2')]}))))

            if tuple(curr_mut) not in seen:  # if this is a unique mutation that has not yet been seen
                master_handle.write('\t'.join(ovals) + '\n')
                seen.add(tuple(curr_mut))
            else:
                sys.stderr.write('Duplicate mutation observed in ' + maf_file + ':\n' + maf_line + '\n')
        maf_handle.close()
    master_handle.close()
    return output_maf


########################################################################################################
# MUTATIONAL BIASES
########################################################################################################

def mutation_per_patient(in_maf_file, cancer_type):
    """
    :param in_maf_file: full path to a .maf-formatted input file with mutations (possibly across multiple samples)
    :param cancer_type: name of the cancer type corresponding to the input maf file
    :return: 2 dictionaries: (1) sample_id -> [# A/T mutations, # C/G mutations] and (2) cancer_type -> set(sample_ids)
    """

    sample_ids = {}  # sample ID -> set(tuples (chromosome, chromosome_position, ref_nucleotide, alt_nucleotide) )
    samples_per_cancer = {cancer_type: set()}  # cancer type -> set( sample IDs )

    mutation_handle = gzip.open(in_maf_file) if in_maf_file.endswith('gz') else open(in_maf_file)

    progress_bar = 1
    header = None
    for mutation_line in mutation_handle:
        if mutation_line.startswith('#'):
            continue
        if not header:
            header = mutation_line.lower()[:-1].split('\t')
            continue

        if progress_bar % 500000 == 0:
            sys.stderr.write('Processed ' + "{:,}".format(progress_bar) + ' lines\n')

        m = mutation_line[:-1].split('\t')
        if m[header.index('start_position')] != m[header.index('end_position')]:  # restrict to SNPs
            continue
        tumor_sample_id = '-'.join(m[header.index('tumor_sample_barcode')].split('-')[:4])
        chromosome = m[header.index('chromosome')]
        chromosome_position = m[header.index('start_position')]
        ref_nucleotide = m[header.index('reference_allele')]
        alt_nuc1 = m[header.index('tumor_seq_allele1')]
        alt_nuc2 = m[header.index('tumor_seq_allele2')]
        alt_nucs = set([v for v in [alt_nuc1, alt_nuc2] if v != ref_nucleotide])

        # keep track of all point mutations for this sample ID
        if tumor_sample_id not in sample_ids:
            sample_ids[tumor_sample_id] = set()
        for alt_nucleotide in alt_nucs:
            sample_ids[tumor_sample_id].add((chromosome, chromosome_position, ref_nucleotide, alt_nucleotide))

        samples_per_cancer[cancer_type].add(tumor_sample_id)  # add this tumor sample to the cancer type

        progress_bar += 1
    mutation_handle.close()

    # we needed to keep all mutations to remove duplicates, now we can summarize:
    biases_per_sample = {}

    for sample_id in sample_ids.keys():

        mutation_counts = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

        for _, _, ref_nuc, alt_nuc in sample_ids[sample_id]:
            if ref_nuc != alt_nuc:  # just a sanity conditional
                mutation_counts[ref_nuc] += 1
        biases_per_sample[sample_id] = [mutation_counts['A'] + mutation_counts['T'],
                                        mutation_counts['C'] + mutation_counts['G']]

    return biases_per_sample, samples_per_cancer


########################################################################################################

def summarize_mutations_by_cancer(in_maffiles, out_biasfile):
    """
    :param in_maffiles: dictionary of cancer_type -> [full path to .maf file, ...]
    :param out_biasfile: full path to write a tab-delimited file of mutational biases to
    :return: none
    """

    biases_per_sample = {}  # sample_id -> [#A/T mutations, #C/G mutations]
    samples_per_cancer = {}  # cancer_type -> set(sample_ids)
    maf_files_per_cancer = {}  # cancer_type -> ','-delimited list of maf files

    # (1) get all the A/T and C/G biases per sample
    for cancer_type, c_maffiles in in_maffiles.items():

        maf_files_per_cancer[cancer_type] = ','.join(c_maffiles)

        for in_maf in c_maffiles:
            current_sample_bias, current_cancer_samples = mutation_per_patient(in_maf, cancer_type)

            # update the current samples
            for sample_id, (at_count, cg_count) in current_sample_bias.items():
                if sample_id not in biases_per_sample:
                    biases_per_sample[sample_id] = [at_count, cg_count]
                else:
                    biases_per_sample[sample_id][0] += at_count
                    biases_per_sample[sample_id][1] += cg_count

            # update the lists of samples per cancer
            for c_type, sample_set in current_cancer_samples.items():
                if c_type not in samples_per_cancer:
                    samples_per_cancer[c_type] = sample_set
                else:
                    samples_per_cancer[c_type] = samples_per_cancer[c_type].union(sample_set)

    # (2) calculate the per-cancer biases:
    cancer_to_bias = {}

    for cancer_type in samples_per_cancer.keys():

        average_ratios = []

        for sample_id in samples_per_cancer[cancer_type]:
            # proportion of overall mutations that were A/T
            sample_ratio = biases_per_sample[sample_id][0] / float(sum(biases_per_sample[sample_id]))
            average_ratios.append(sample_ratio)

        average_ratios = sum(average_ratios) / len(average_ratios)

        cancer_to_bias[cancer_type] = average_ratios

    # (2) for each cancer, print out the total bias
    bias_handle = gzip.open(out_biasfile, 'wt') if out_biasfile.endswith('gz') else open(out_biasfile, 'w')
    bias_handle.write('\n'.join(['# Mutational biases evaluated from single nucleotide point mutations',
                                 '# Rate of A/T vs. C/G mutations were assessed *per tumor sample* and averaged ' +
                                 'by cancer type',
                                 '\t'.join(['#cancer_type', 'A/T_mutation_frequency', 'C/G_mutation_frequency',
                                            'N_mutation_frequency', 'total_samples', 'maf_file'])]) + '\n')

    for average_ratio, cancer_type in sorted([(ratio, cancer_type) for cancer_type, ratio in cancer_to_bias.items()]):
        cg_likelihood = (1 - average_ratio) / average_ratio
        bias_handle.write('\t'.join([cancer_type,
                                     '1',
                                     str(cg_likelihood),
                                     str((cg_likelihood + 1) / 2.),
                                     str(len(samples_per_cancer[cancer_type])),
                                     maf_files_per_cancer[cancer_type]]) + '\n')
        sys.stderr.write(cancer_type + '\t' + "{:.2f}".format(cg_likelihood) + 'x more likely C/G mutation\n')
    bias_handle.close()

    sys.stderr.write('Mutational biases written to ' + out_biasfile + '\n')


########################################################################################################

def process_mutational_bias(mutation_bias_file, cancer_type):
    """
    :param mutation_bias_file: full path to a tab-delimited file (produced by mutation_biases.py) with columns
                               corresponding to cancer name, relative likelihood of an A/T mutation, relative
                               likelihood of a C/G mutation and relative likelihood of an N mutation (average of
                               A/T and C/G)
    :param cancer_type: string corresponding to a specific cancer type to get mutation biases for
    :return: dictionary of nucleotide -> relative mutational likelihood
    """

    if not os.path.isfile(mutation_bias_file):
        sys.stderr.write('Could not open ' + mutation_bias_file + '\n')
        sys.exit(1)

    biases = {'A': 1., 'a': 1., 'T': 1., 't': 1.,
              'C': 1., 'c': 1., 'G': 1., 'g': 1.,
              'N': 1., 'n': 1.}

    bias_handle = gzip.open(mutation_bias_file) if mutation_bias_file.endswith('gz') else open(mutation_bias_file)
    for bias_line in bias_handle:
        if bias_line.startswith('#'):
            continue
        current_cancer_type = bias_line[:-1].split('\t')[0]
        at_mutation, cg_mutation, n_mutation = map(float, bias_line[:-1].split('\t')[1:4])

        if current_cancer_type == cancer_type:
            for nucleotide_set, observed_bias in [(('A', 'a', 'T', 't'), at_mutation),
                                                  (('C', 'c', 'G', 'g'), cg_mutation),
                                                  (('N', 'n'), n_mutation)]:
                for nucleotide in nucleotide_set:
                    biases[nucleotide] = observed_bias
            input_files = bias_line[:-1].split('\t')[5]
            break
    else:
        sys.stderr.write('Could not find mutational biases for ' + cancer_type + ' in ' + mutation_bias_file + '\n')
        sys.stderr.write('Returning equal likelihood of mutation at each nucleotide basepair...\n')

    return biases, input_files


########################################################################################################

def missense_likelihood(codon, codon_position):
    """
    :param codon: 3-letter string composed of As, Ts, Cs, Gs, Ns
    :param codon_position: index (0, 1, 2) of a CHANGE at a codon position
    :return: the likelihood (either 0/3, 1/3, 2/3, or 3/3) of a mutation at the specified codon position
             to cause a missense mutation
    """

    # hardcode the mapping from codon -> amino acid
    codon_to_amino_acid = GENCODE

    start_codon = codon.upper()
    if start_codon not in codon_to_amino_acid:  # unknown how nucleotide changes cause AA changes
        return 1.

    start_amino_acid = codon_to_amino_acid[start_codon]
    change_count = 0  # for all 4 possible nucleotide changes, how many resulted in a missense mutation?
    for changed_nucleotide in ['A', 'T', 'C', 'G']:
        if changed_nucleotide != start_codon[codon_position]:
            new_codon = ''.join([current_nucleotide if codon_index != codon_position else changed_nucleotide for
                                 codon_index, current_nucleotide in enumerate(list(start_codon))])
            if new_codon not in codon_to_amino_acid or codon_to_amino_acid[new_codon] != start_amino_acid:
                change_count += 1
    return change_count / (3. if start_codon[codon_position] != 'N' else 4.)


########################################################################################################

def get_mutation_likelihood(cancer_type, cdna_file, mutation_bias_file, final_frequency_file, bias_type='both'):
    """
    :param cancer_type: type of cancer to obtain C/G mutation biases for
    :param cdna_file: full path to a fasta-formatted file with cDNA sequences (to use to create per-position
                      mutation likelihood "tracks" per protein)
    :param mutation_bias_file: full path to the file created by summarize_mutations_by_cancer()
    :param final_frequency_file: full path to an output file to write results to
    :param bias_type: either "cg", "missense", or "both" for how to evaluate the likelihood of each protein position
                      to harbor a missense mutation
    :return: none
    """

    # go through the cDNA file, keeping track of the sequence ID and also the likelihood of mutation (not relative)
    #  at each amino acid position:
    if not os.path.isfile(cdna_file):
        sys.stderr.write('Could not open ' + cdna_file + '\n')
        sys.exit(1)

    if bias_type in ['gc', 'both']:
        biases, input_maf = process_mutational_bias(mutation_bias_file, cancer_type)

    mutation_frequency_file = final_frequency_file.replace('.txt', '_' + cancer_type + '.txt')
    if mutation_frequency_file.endswith('gz'):
        freq_outhandle = gzip.open(mutation_frequency_file, 'w')
    else:
        freq_outhandle = open(mutation_frequency_file, 'w')

    freq_outhandle.write('\n'.join(['# Mutational biases *for ' + cancer_type + ' cancer* evaluated from single ' +
                                    'nucleotide point mutations found in',
                                    '# ' + input_maf,
                                    '# All cDNA sequences found in ' + cdna_file,
                                    '\t'.join(['#sequence_id', 'protein_length', 'AN(irrelevant)',
                                               '0-index-position:mutation_likelihood,...'])]) + '\n')

    cdna_handle = gzip.open(cdna_file) if cdna_file.endswith('gz') else open(cdna_file)
    for fasta_line in cdna_handle:
        if fasta_line.startswith('>'):
            seq_id = fasta_line[1:-1].split()[0]
            cdna_seq = cdna_handle.next().strip()

            likelihoods = []
            for codon_start_index in xrange(0, len(cdna_seq), 3):
                codon = cdna_seq[codon_start_index:codon_start_index + 3]  # current codon
                if len(codon) < 3:
                    continue

                if bias_type in ['gc', 'both']:
                    mutation_likelihoods = [biases[codon[codon_position]] for codon_position in xrange(3)]
                else:
                    mutation_likelihoods = [1., 1., 1.]  # the chance of each nucleotide being hit (GC bias)

                if bias_type in ['missense', 'both']:
                    missense_likelihoods = [missense_likelihood(codon, codon_position) for codon_position in xrange(3)]
                else:
                    missense_likelihoods = [1., 1., 1.]  # the chance of each nucleotide causing a missense mutation

                # the overall likelihood of a missense mutation is therefore:
                current_missense_likelihood = sum([mutation_likelihoods[i] * missense_likelihoods[i]
                                                   for i in xrange(3)])
                likelihoods.append(current_missense_likelihood)

            # finished processing this sequence, write to file:
            # NOTE: these lambdas *do not need to sum to 1* -- they are normalized with respect to each track when
            #       calculating the expectation
            freq_outhandle.write('\t'.join([seq_id,  # sequence ID
                                            str(int(len(cdna_seq) / 3.)),  # length of protein
                                            '0',  # irrelevant (vestigial)
                                            ','.join([str(prot_index) + ':' + str(mut_likelihood) for
                                                      prot_index, mut_likelihood in enumerate(likelihoods)])]) + '\n')
    cdna_handle.close()
    freq_outhandle.close()
    sys.stderr.write('Finished processing mutational biases! File in ' + mutation_frequency_file + '\n')


########################################################################################################
# PROCESS INPUT DATA
########################################################################################################

def get_maf_file(mutation_path, cancer_type, exit_on_error=False):
    """
    :param mutation_path: path to a directory with subdirectories specified by cancer type
    :param cancer_type: type of cancer we are looking at
    :param exit_on_error: boolean indicating whether to exit if a maf file was not found or return without error
    :return: full path to the .maf file corresponding to this cancer type, if possible
    """

    # get the .maf file that we are interested in:
    maf_file = None
    if cancer_type in os.listdir(mutation_path):
        mutation_file = [mutation_path + cancer_type + '/' + file_name
                         for file_name in os.listdir(mutation_path + cancer_type)
                         if file_name.startswith('TCGA.' + cancer_type + '.muse.')
                         and '.somatic.maf' in file_name]
        if len(mutation_file) > 0:
            maf_file = mutation_file[0]

    if not maf_file:
        if exit_on_error:
            sys.stderr.write('No input .maf file found in ' + mutation_path + cancer_type + '/\n' +
                             'Please run: python ' + sys.argv[0] + ' --rename_files\n' +
                             ('python ' + sys.argv[0] + ' --aggregate_files\n'
                              if cancer_type in {'Aggregate', 'AggregateTrunc', 'PANGI', 'KIPAN', 'COADREAD',
                                                 'GBMLGG', 'STES'} else ''))
            sys.exit(1)

    return maf_file


########################################################################################################

def extract_mutations_from_maf(in_maffile, out_mutfile):
    """
    :param in_maffile: input .maf file (http://software.broadinstitute.org/software/igv/MutationAnnotationFormat)
    :param out_mutfile: full path to an output file to write all within-exon mutations to (for further analysis)
    :return: none
    """

    if not os.path.isfile(in_maffile):
        sys.stderr.write('Could not open ' + in_maffile + '...Exiting.\n')
        sys.exit(1)

    all_mutations = set()  # save the mutations that we are seeing
    header = None

    maf_handle = gzip.open(in_maffile) if in_maffile.endswith('.gz') else open(in_maffile)
    for l in maf_handle:
        if l.startswith('#'):
            continue
        if not header:  # lowercase header due to version differences in .maf files
            header = l[:-1].lower().split('\t') if '\t' in l else l[:-1].lower().split()
            continue

        # all the values for the current line:
        m = l[:-1].split('\t') if '\t' in l else l[:-1].split()

        # check required fields:
        if len(m) < max([header.index(i) for i in
                         ['chromosome', 'start_position', 'end_position', 'variant_type', 'variant_classification',
                          'reference_allele',
                          'tumor_seq_allele1', 'tumor_seq_allele2', 'hugo_symbol', 'tumor_sample_barcode']]) + 1:
            continue

        if 'codon_change' not in header and 'hgvsc' not in header:
            sys.stderr.write('Cannot find codon change!\n')
            continue

        if 'protein_change' not in header and 'hgvsp_short' not in header:
            sys.stderr.write('Cannot find protein change!\n')
            continue

        # Extract EXONIC variants only:
        # Non-exon variant classifications =
        # {"3'UTR", "3'Flank", "5'UTR", "5'Flank", 'Intron', 'RNA', 'lincRNA', 'Targeted_Region', 'IGR'}
        if m[header.index('variant_classification')] not in \
            ['De_novo_Start_InFrame', 'De_novo_Start_OutOfFrame', 'Frame_Shift_Del', 'Frame_Shift_Ins', 'In_Frame_Del',
             'In_Frame_Ins', 'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation', 'Silent', 'Splice_Site',
             'Start_Codon_Del', 'Start_Codon_Ins', 'Start_Codon_SNP', 'Stop_Codon_Del', 'Stop_Codon_Ins',
             'Translation_Start_Site']:
            continue

        # extract all alternate identifiers for this gene/protein
        ids = {}
        for seq_id, header_id in [('hugo_symbol', ('hugo_symbol', '')),
                                  ('entrez', ('entrez_gene_id', '')),
                                  ('refseq', ('refseq_prot_id', 'refseq')),  # refseq
                                  ('ensembl_transcript', ('annotation_transcript', 'transcript_id')),  # transcript_id
                                  ('gene_id', ('i_hgnc_ensembl_gene_id', 'gene')),  # gene
                                  ('ensembl_prot', ('i_hgvs_protein_change', 'ensp'))]:  # ensp
            idmatches = [m[header.index(hid)] for hid in header_id if hid in header and len(hid) > 0]
            ids[seq_id] = idmatches[0] if len(idmatches) > 0 else ''
        ids['ensembl_prot'] = ids['ensembl_prot'].split(':')[0]

        gene_id = ids['gene_id']  # change to hugoSymbol, entrez, refseq, etc. as desired

        # Extract location, cDNA and amino acid changes:
        genome_build = m[header.index('ncbi_build')]
        chrom_loc = m[header.index('chromosome')]
        if chrom_loc.startswith('chr'):
            chrom_loc = chrom_loc[3:]
        start_pos = m[header.index('start_position')]
        end_pos = m[header.index('end_position')]
        mut_type = m[header.index('variant_classification')].replace('_Mutation', '')
        ref_allele = m[header.index('reference_allele')]
        alt_allele = m[header.index('tumor_seq_allele1')] + '|' + m[header.index('tumor_seq_allele2')]

        if 'codon_change' in header:
            codon_index_val = 'codon_change'
        elif 'hgvsc' in header:
            codon_index_val = 'hgvsc'

        curr_codon = m[header.index(codon_index_val)][m[header.index(codon_index_val)].find(')') + 1:]

        if '>' in curr_codon and (codon_index_val == 'codon_change' or 'codons' in header):
            ref_codon = curr_codon.split('>')[0] if codon_index_val == 'codon_change' else \
                m[header.index('codons')].split('/')[0]
            try:
                alt_codon = curr_codon.split('>')[1] if codon_index_val == 'codon_change' else \
                    m[header.index('codons')].split('/')[1]
            except:
                alt_codon = ''
            codon_pos = ','.join([str(i + 1) for i, a in enumerate(list(ref_codon)) if a.isupper()])
        elif 'fs' in curr_codon:
            ref_codon = curr_codon.replace('fs', '')
            alt_codon = curr_codon.replace('fs', '')
            codon_pos = 'fs'
        else:
            ref_codon = m[header.index(codon_index_val)][m[header.index(codon_index_val)].find('c.') + 2:]
            alt_codon = m[header.index(codon_index_val)][m[header.index(codon_index_val)].find('c.') + 2:]
            codon_pos = m[header.index(codon_index_val)][m[header.index(codon_index_val)].find('c.') + 2:]

        if 'protein_change' in header:
            protein_index_val = 'protein_change'
        elif 'hgvsp_short' in header:
            protein_index_val = 'hgvsp_short'
        curr_prot = m[header.index(protein_index_val)][m[header.index(protein_index_val)].find('p.') + 2:]

        if '>' in curr_prot:
            prot_pos = ''.join([a for a in curr_prot.split('(')[0].split('_')[0] if a.isdigit()])
            ref_prot = ''.join([a for a in curr_prot.replace('_', '').replace('(', '').replace(')', '')
                                if not a.isdigit()]).split('>')[0]
            alt_prot = ''.join([a for a in curr_prot.replace('_', '').replace('(', '').replace(')', '')
                                if not a.isdigit()]).split('>')[1]
        elif len(curr_prot) > 0:
            numbers = [i for i, a in enumerate(list(curr_prot)) if is_numeric(a)]
            if len(numbers) > 0:
                prot_pos = curr_prot[min(numbers):max(numbers) + 1]
                ref_prot = curr_prot[:min(numbers)]
                alt_prot = curr_prot[max(numbers) + 1:]
            else:
                prot_pos, ref_prot, alt_prot = '', '', ''
        else:
            prot_pos, ref_prot, alt_prot = '', '', ''

        sample_id = '-'.join(m[header.index('tumor_sample_barcode')].split('-')[:4])

        all_mutations.add((gene_id,
                           ','.join([seq_id + ':' + ids[seq_id] for seq_id in
                                     ['hugo_symbol', 'entrez', 'refseq', 'gene_id', 'ensembl_transcript',
                                      'ensembl_prot'] if seq_id in ids and ids[seq_id] != '']),
                           sample_id,
                           mut_type,
                           genome_build,
                           chrom_loc,
                           start_pos + '-' + end_pos if start_pos != end_pos else start_pos,
                           ref_allele,
                           alt_allele,
                           codon_pos,
                           ref_codon,
                           alt_codon,
                           prot_pos,
                           ref_prot,
                           alt_prot))
    maf_handle.close()

    out_handle = gzip.open(out_mutfile, 'w') if out_mutfile.endswith('.gz') else open(out_mutfile, 'w')
    out_handle.write('# All exon mutations extracted from: ' + in_maffile + '\n')
    out_handle.write('\t'.join(
        ['#SequenceID', 'OtherIDs', 'SampleID', 'MutationType', 'NCBI_Build', 'Chromosome', 'ChromPos', 'RefAllele',
         'AltAllele', 'CodonPos', 'RefCodon', 'AltCodon', 'ProtPos', 'RefAA', 'AltAA']) + '\n')
    out_handle.write('\n'.join(['\t'.join(a) for a in sorted(list(all_mutations))]) + '\n')
    out_handle.close()
    sys.stderr.write('Condensed ' + in_maffile + ' into ' + out_mutfile + '\n')


########################################################################################################
# FIND MUTATIONS FALLING INTO PROTEIN-CODING EXONS
########################################################################################################

def find_exon_locations(fasta_file):
    """
    :param fasta_file: Path to directory of Ensembl genes
    :return: an IntervalTree of gene locations on the appropriate chromosome
    """

    if not os.path.isfile(fasta_file):
        sys.stderr.write('Could not open ' + fasta_file + '\n')
        return {}

    gene_locations = {}
    gene_intervals = {}

    # extract exon start/end information for all exons by gene
    seq_handle = gzip.open(fasta_file) if fasta_file.endswith('gz') else open(fasta_file)
    for seq_line in seq_handle:
        if seq_line.startswith('>'):

            loc = seq_line[seq_line.find(BUILD + ':') + len(BUILD + ':'):].split()[0]
            chrom_id = loc.split(':')[0]
            if chrom_id.startswith('chr'):
                chrom_id = chrom_id[3:]

            exons = loc.split(':')[1].replace('join(', '').replace('complement(', '').replace(')', '').split(',')
            exons = [tuple(sorted([int(a.split('..')[0]), int(a.split('..')[1])])) for a in exons]

            gene_id = seq_line[seq_line.find('gene:') + 5:].split()[0]  # e.g., ENSG00000000457.9

            if chrom_id not in gene_intervals:
                gene_intervals[chrom_id] = {}
            if gene_id not in gene_intervals[chrom_id]:
                gene_intervals[chrom_id][gene_id] = set()

            # Add each POSITIVE exon SEPARATELY!
            for exon_start, exon_end in exons:  # these are 1-indexed, but we need 0-indexed...
                if exon_start > -1 and exon_end > -1:
                    gene_intervals[chrom_id][gene_id].add((exon_start, exon_end))
    seq_handle.close()

    # Store all the information as interval trees:
    for chrom_id in gene_intervals.keys():
        for gene_id in gene_intervals[chrom_id].keys():
            gene_min_start = sys.maxint
            gene_max_end = -1

            for exon_start, exon_end in gene_intervals[chrom_id][gene_id]:
                if exon_start < gene_min_start:
                    gene_min_start = exon_start
                if exon_end > gene_max_end:
                    gene_max_end = exon_end

            if chrom_id not in gene_locations:
                gene_locations[chrom_id] = IntervalTree()

            for exon_start, exon_end in gene_intervals[chrom_id][gene_id]:
                gene_locations[chrom_id].addi(exon_start,
                                              exon_end + (1 if exon_end == exon_start else 0),
                                              (gene_id, gene_min_start, gene_max_end))

    return gene_locations


########################################################################################################

def find_correct_gene_from_location(gene_id, chrom, origchrompos, genelocs):
    """Given a genomic location, determine the gene, if any, that the mutation
    falls into (provided that the original gene was incorrect)."""

    if chrom not in genelocs:
        return ''

    all_matching_genes = set()

    chrompos_start = int(origchrompos.split('-')[0])
    chrompos_end = int(origchrompos.split('-')[-1]) + 1
    for start, end, exon_id in genelocs[chrom][chrompos_start:chrompos_end]:
        gene = exon_id[0]

        if gene != gene_id:
            all_matching_genes.add(gene)

    if len(all_matching_genes) == 1:
        return list(all_matching_genes)[0]
    elif len(all_matching_genes) == 0:
        return ''
    else:
        sys.stderr.write('multiple matching genes: ' + ', '.join(list(all_matching_genes)) + '\n')
        return list(all_matching_genes)[0]


########################################################################################################

def get_alternate_ids(header, other_ids):
    """
    :param header: sequence header with alternate IDs
    :param other_ids: set of alternate IDs that have already been extracted
    :return: update alternate IDs from some input file based on what we found from Ensembl's BioMart
    """

    ensembl_gene_id = header[header.find('gene:') + 5:].strip().split()[0] if 'gene' in header else ''
    ensembl_cdna_id = header[header.find('transcript:') + 11:].strip().split()[0] if 'transcript' in header else ''
    ensembl_prot_id = header[header.find('prot:') + 5:].strip().split()[0] if 'prot' in header else ''
    entrez_gene_id = header[header.find('entrez:') + 7:].strip().split()[0] if 'entrez' in header else \
        (other_ids[other_ids.find('entrez:') + 7:].strip().split(',')[0] if 'entrez' in other_ids else '')
    refseq_prot_id = header[header.find('refseq:') + 7:].strip().split()[0] if 'refseq' in header else \
        (other_ids[other_ids.find('refseq:') + 7:].strip().split(',')[0] if 'refseq' in other_ids else '')
    hgnc_id = header[header.find('hgncID:') + 7:].strip().split()[0] if 'hgncID' in header else ''
    hugo_gene_id = header[header.find('hugoSymbol:') + 11:].strip().split()[0] if 'hugoSymbol' in header else \
        (other_ids[other_ids.find('hugoSymbol:') + 11:].strip().split(',')[0] if 'hugoSymbol' in other_ids else '')

    return {'gene_id': ensembl_gene_id,
            'ensembl_transcript': ensembl_cdna_id,
            'ensembl_prot': ensembl_prot_id,
            'entrez': entrez_gene_id,
            'refseq': refseq_prot_id,
            'hugoSymbol': hugo_gene_id,
            'hgncID': hgnc_id}


########################################################################################################
# VERIFY AND FIX MUTATIONS
########################################################################################################

def find_frameshifts(mut_infile):
    """Return a set of frameshift mutation locations, per sample, so that
    downstream mutations are ignored in our analysis."""

    sys.stderr.write('Finding frameshift mutations in ' + mut_infile + '...\n')
    if not os.path.isfile(mut_infile):
        sys.stderr.write('Could not open ' + mut_infile + ', returning empty list of frameshift mutations.\n')
        return {}

    frame_shifts = {}  # sampleID -> set of tuples: (chromosome, chromosome position, mutation type)
    exon_inhandle = gzip.open(mut_infile) if mut_infile.endswith('.gz') else open(mut_infile)
    for exon_mutline in exon_inhandle:
        if exon_mutline.startswith('#') or len(exon_mutline[:-1].split('\t')) < 7:
            continue
        seq_id, other_ids, sample_id, mut_type, genome_build, chrom, pos = exon_mutline.strip().split('\t')[:7]
        if mut_type not in ['Missense', 'Silent']:
            if sample_id not in frame_shifts:
                frame_shifts[sample_id] = set()
            frame_shifts[sample_id].add((chrom, pos, mut_type))
    exon_inhandle.close()

    return frame_shifts


########################################################################################################

def edit_mutations_by_line(l, genelocs, frameshifts, outfile_handle, logfile_handle,
                           exon_dir=data_path + 'ensembl/Homo_sapiens.' + BUILD + 'exons/',
                           remove_frameshifts=True):
    """
    :param l: line from an exonic mutation file
    :param genelocs: exon intervals mapping to gene identifiers
    :param frameshifts: sampleID -> [(chrom, chrom_pos)] of frameshift mutations
    :param outfile_handle: output file handle (already opened for writing) to write newly-corrected lines
    :param logfile_handle: output file handle (already opened for writing) to write errors
    :param exon_dir: directory to all exon cDNA information for each protein isoform
    :param remove_frameshifts: boolean indicating whether we want to remove mutations after frameshifts or not
    :return: none
    """

    if len(l[:-1].split('\t')) < 15:
        logfile_handle.write(l[:-1] + '\tERROR:TOO_SHORT\n')
        return

    (gene_id, other_ids, sample_id, mut_type, genome_build, chrom, orig_chrom_pos, orig_ref_allele, orig_alt_allele,
     codon_pos, ref_codon, alt_codon, prot_pos, ref_prot, alt_prot) = l[:-1].split('\t')[:15]

    if chrom.startswith('chr'):
        chrom = chrom[3:]

    # ['37','GRCh37','GRCh37-lite','hg19']: # '36', '36.1', 'hg18'
    if genome_build not in [BUILD, BUILD + '-lite', BUILD.replace('GRCh', ''), 'hg' + BUILD.replace('GRCh', '')]:
        logfile_handle.write(l[:-1] + '\tERROR:INCORRECT_BUILD__' + genome_build + '\n')
        return

    # 'Start_Codon_SNP', 'Start_Codon_Del', 'Splice_Site', 'De_novo_Start_InFrame'
    if mut_type not in ['Missense', 'Silent', 'Nonsense', 'Nonstop']:
        logfile_handle.write(l[:-1] + '\tERROR:OTHER_MUTATION_TYPE\n')
        if remove_frameshifts:
            return

    if int(orig_chrom_pos.split('-')[-1]) - int(orig_chrom_pos.split('-')[0]) + 1 != len(
        orig_ref_allele.replace('-', '')):
        logfile_handle.write(l[:-1] + '\tERROR:CHROMOSOME_POSITION_INCORRECT_' + \
                             str(int(orig_chrom_pos.split('-')[-1]) - int(orig_chrom_pos.split('-')[0]) + 1) + '!=' + \
                             str(len(orig_ref_allele.replace('-', ''))) + '\n')
        return

    # Get the complete Ensembl gene ID if necessary (our data is 1-to-1)
    if gene_id.startswith('ENSG') and '.' not in gene_id and chrom in os.listdir(exon_dir):
        # e.g., ENSG00000132698 -> ENSG00000132698.9
        gene_ids = sorted([a for a in os.listdir(exon_dir + chrom) if a.startswith(gene_id)])
        if len(gene_ids) > 0:
            gene_id = gene_ids[0]  # only relevant for some ambiguous user input

    if chrom == 'M':
        chrom = 'MT'  # Relabel any mitochrondial sequences to match exon labels

    #### In some cases, the gene ID specified in mutation file does not match the
    #### location we have for corresponding Ensembl genes, so we attempt to edit:
    if (chrom not in os.listdir(exon_dir) or
        gene_id == '' or
        gene_id not in os.listdir(exon_dir + chrom) or
        len(os.listdir(exon_dir + chrom + '/' + gene_id)) < 1):

        new_gene_id = find_correct_gene_from_location(gene_id, chrom, orig_chrom_pos, genelocs)

        if (new_gene_id in [gene_id, ''] or
            chrom not in os.listdir(exon_dir) or
            new_gene_id not in os.listdir(exon_dir + chrom) or
            len(os.listdir(exon_dir + chrom + '/' + new_gene_id)) < 1):
            logfile_handle.write(l[:-1] + '\tERROR:IRRETRIEVABLE_GENE_ID\n')
            return
        else:
            gene_id = new_gene_id

    # Determine if there are any possible frameshift mutations in this sample, in this gene
    frameshift_changes = set()  # all positions *within this gene* corresponding to frameshifts/truncatations/etc.
    if sample_id in frameshifts:  # frameshifts -> sampleID -> [(chrom, chrom_pos)]
        for fs_chrom, fs_chrompos, fs_muttype in frameshifts[sample_id]:
            if fs_chrom == chrom:
                frameshift_start = int(fs_chrompos.split('-')[0])
                frameshift_end = int(fs_chrompos.split('-')[-1]) + 1
                for frameshift_loc in range(frameshift_start, frameshift_end):
                    for exon_start, exon_end, exon_id in genelocs[chrom][frameshift_loc]:
                        this_gene = exon_id[0]
                        if this_gene == gene_id:
                            frameshift_changes.add(frameshift_loc)

    # Consider EACH SNP (e.g., also from DNPs, TNPs, ONPs) individually:
    if '|' not in orig_alt_allele:
        orig_alt_allele = orig_alt_allele + '|' + orig_alt_allele
    # if '|' not in orig_alt_allele: orig_alt_allele = orig_ref_allele+'|'+orig_alt_allele

    for k in xrange(max(len(orig_ref_allele),
                        len(orig_alt_allele.split('|')[0]),
                        len(orig_alt_allele.split('|')[1]))):
        # chrom_pos = range(int(orig_chrom_pos.split('-')[0]), int(orig_chrom_pos.split('-')[-1]) + 1)[k]
        chrom_pos = int(orig_chrom_pos.split('-')[0]) + k
        ref_allele = orig_ref_allele[k] if len(orig_ref_allele) > k else '-'
        alt_allele = (orig_alt_allele.split('|')[0][k] if len(orig_alt_allele.split('|')[0]) > k else '-') + '|' + \
                     (orig_alt_allele.split('|')[1][k] if len(orig_alt_allele.split('|')[1]) > k else '-')

        if ref_allele == '-' or alt_allele == '-|-':
            logfile_handle.write(l[:-1] + '\tERROR:NOT_A_SNP\n')
            return

        #### We match each SNP to the reference genome, then incorporate the change and
        #### translate to functionally classify the SNP
        for exon_file in sorted([a for a in os.listdir(exon_dir + chrom + '/' + gene_id) if a.endswith('.exons.txt')]):
            curr_prot_id = exon_file.replace('.exons.txt', '')  # e.g., ENSP00000476288.1
            exons = []
            with open(exon_dir + chrom + '/' + gene_id + '/' + exon_file) as y:
                alt_ids = get_alternate_ids(y.next(), other_ids)
                for exon in y:
                    pos = map(int, exon.strip().split('\t')[0].split(':'))  # start -> end of exon
                    seq = exon.strip().split('\t')[1]  # sequence, reversed + strand in case of complement
                    exons.append((tuple(pos), seq))

            # --------------------------------------------------------------------------------------------------------------
            #### If there ARE frameshifts, keep track of the [0] exon rank, [1] position in exon, and
            #### [2] genomic location of frameshift mutations
            problems = set()
            if len(frameshift_changes) > 0:
                for exon_rank, (pos, seq) in enumerate(exons):
                    complement = pos[0] > pos[1]
                    this_exon_pos = range(min(pos), max(pos) + 1)[::-1 if complement else 1]
                    for frameshift_loc in [(exon_rank, this_exon_pos.index(a), a) for a in frameshift_changes if
                                           a in this_exon_pos]:
                        problems.add(frameshift_loc)
            problems = sorted(list(problems))  # ordered by exon rank, exon location, genomic location

            # Finally, go through each exon to determine where 'chrom_pos' should be:
            for exon_rank, (pos, seq) in enumerate(exons):
                complement = pos[0] > pos[1]
                this_exon_pos = range(min(pos), max(pos) + 1)[::-1 if complement else 1]

                if chrom_pos in this_exon_pos:  # THIS is the exon where our mutation has landed!
                    this_index = this_exon_pos.index(chrom_pos)
                    if ref_allele == 'N':
                        ref_allele = seq[this_index].upper()
                    if not seq[this_index].upper() == ref_allele.upper() or ref_allele == '-':
                        logfile_handle.write(l[:-1] + '\tERROR:DID_NOT_MATCH_REFERENCE_' + ref_allele.upper() + '->' +
                                             seq[this_index].upper() + '\n')
                        break  # try the next protein sequence

                    if len(problems) > 0 and (exon_rank > problems[0][0] or
                                              (exon_rank == problems[0][0] and this_index >= problems[0][1])):
                        logfile_handle.write(l[:-1] + '\tERROR:AFTER_A_FRAMESHIFT\n')
                        if remove_frameshifts:
                            break

                    # Extract codon information
                    prev_cdna_len = sum(
                        [len(exons[j][1]) for j in xrange(exon_rank)])  # other exons in order up to this point
                    this_cdna_index = prev_cdna_len + this_index  # and the length within this exon..
                    this_codon_index = ((this_cdna_index + 1) / 3) + (1 if (this_cdna_index + 1) % 3 > 0 else 0)
                    pos_in_codon = (this_cdna_index + 1) % 3 if (this_cdna_index + 1) % 3 > 0 else 3

                    try:
                        if pos_in_codon == 1:
                            this_codon = seq[this_index:this_index + 3]
                        elif pos_in_codon == 2:
                            this_codon = seq[this_index - 1:this_index + 2]
                        elif pos_in_codon == 3:
                            this_codon = seq[this_index - 2:this_index + 1]
                        if len(this_codon) != 3: raise IndexError
                    except IndexError:  # need to get into the previous or next exon to get the actual codon..
                        new_seq, new_this_index = seq, this_index
                        if exon_rank > 0:  # not the first exon, so look to the left
                            new_seq = exons[exon_rank - 1][1] + new_seq
                            new_this_index = this_index + len(exons[exon_rank - 1][1])
                        if exon_rank < len(exons) - 1:  # not the last exon, so also look to the right
                            new_seq = new_seq + exons[exon_rank + 1][1]
                        if pos_in_codon == 1:
                            this_codon = new_seq[new_this_index:new_this_index + 3]
                        elif pos_in_codon == 2:
                            this_codon = new_seq[new_this_index - 1:new_this_index + 2]
                        elif pos_in_codon == 3:
                            this_codon = new_seq[new_this_index - 2:new_this_index + 1]

                        if len(this_codon) != 3:
                            logfile_handle.write(l[:-1] + '\tERROR:UNABLE_TO_COMPLETE_CODON\n')
                            break

                    # Get the new reference codon (capitalize mutation position) and newly translated reference protein
                    nref_codon = ''.join([COMPLEMENT[a] if complement else a for a in this_codon]).lower()
                    nref_codon = ''.join([a.upper() if j == pos_in_codon - 1 else a for j, a in enumerate(nref_codon)])
                    nref_prot = GENCODE.get(nref_codon.upper(), 'X')

                    # And now compare that to each alternate codon and translated amino acid
                    for d in sorted(list(set([d for d in alt_allele.split('|') if d != ref_allele and d != '-']))):
                        nalt_codon = ''.join([nref_codon[j] if j != pos_in_codon - 1 else
                                              (COMPLEMENT[d] if (complement and 'GWAS' not in sample_id) else d) for j
                                              in
                                              xrange(3)])
                        if nref_codon.lower() == nalt_codon.lower():
                            continue

                        nalt_prot = GENCODE.get(nalt_codon.upper(), 'X')
                        nmut_type = 'Silent' if nref_prot == nalt_prot else (
                            'Nonstop' if nref_prot == '_' else ('Nonsense' if nalt_prot == '_' else 'Missense'))

                        new_mut_line = '\t'.join([curr_prot_id,
                                                  ','.join([i + ':' + alt_ids[i] for i in
                                                            ['hugoSymbol', 'hgncID', 'entrez', 'refseq', 'gene_id',
                                                             'ensembl_transcript', 'ensembl_prot'] if
                                                            alt_ids[i] != '']),
                                                  sample_id,
                                                  nmut_type,
                                                  genome_build,
                                                  chrom,
                                                  str(chrom_pos),
                                                  ref_allele,
                                                  d,
                                                  str(pos_in_codon),
                                                  nref_codon, nalt_codon,
                                                  str(this_codon_index), nref_prot, nalt_prot]) + '\n'

                        if nmut_type != mut_type:  # our assessment of the mutation type is INCORRECT!
                            logfile_handle.write(l[:-1] + '\tERROR:INCORRECT_MUTATION_TYPE_' +
                                                 mut_type + '->' + nmut_type + '\n')
                            # continue # It seems we definitely want to keep these, anyway!

                        outfile_handle.write(new_mut_line)


########################################################################################################

def remap_mutations_to_all_proteins(mut_infile, verified_outfile, log_file='',
                                    seq_file=data_path + 'ensembl/Homo_sapiens.GRCh38/' +
                                             'Homo_sapiens.GRCh38.pep.all.withgenelocs.verified.fa.gz',
                                    exon_dir=data_path + 'ensembl/Homo_sapiens.GRCh38/exons/',
                                    inc_frameshifts=False):
    """
    :param mut_infile: full path to a file generated by extract_mutations_from_maf function
      > [0] sequence_id, [1] alternate_ids, [2] sample_id, [3] mutation_type, [4] genome_build, [5] chromosome,
        [6] chrom_position, [7] ref_allele, [8] alt_allele, [9] codon_offset, [10] ref_codon, [11] alt_codon,
        [12] prot_position, [13] ref_aa, [15] alt_aa
    :param seq_file: full path to a fasta file with location information in the sequence headers (to quickly
                     identify genes to search)
    :param verified_outfile: full path to write verified mutation output to
    :param log_file: full path to write errors to
    :param exon_dir: full path to a directory containing subdirectories separated by chromosome name
                     with .exons.txt files in each Ensembl gene subdirectory
    :param inc_frameshifts: boolean indicating whether to include mutations after frameshift/nonsense mutations
    :return: none, but map mutations onto ALL protein isoforms to verify variant calls (there may be multiple
             isoforms per gene)
    """

    # (1) Read in gene locations in order to remap and edit mutations, line by line
    gene_locations = find_exon_locations(seq_file)

    # (2) Identify all potential FRAMESHIFT mutations, to avoid reporting those mutations
    #     that occur after frameshifts in a per-sample basis
    frameshifts = {}
    if not inc_frameshifts:
        frameshifts = find_frameshifts(mut_infile)  # sampleID -> [(chrom, chrom_pos)]

    # (3) rename error file if left unspecified
    if log_file == '':
        log_file = verified_outfile.replace('.txt', '-ERRORS.txt')

    out_verified = gzip.open(verified_outfile, 'w') if verified_outfile.endswith('.gz') else open(verified_outfile, 'w')
    out_verified.write('# All mutations edited from ' + mut_infile + '\n')
    out_verified.write('# See logfile ' + log_file + ' for all mapping errors\n')
    out_errors = gzip.open(log_file, 'w') if log_file.endswith('.gz') else open(log_file, 'w')

    mut_inhandle = gzip.open(mut_infile) if mut_infile.endswith('.gz') else open(mut_infile)
    for mut_line in mut_inhandle:
        if mut_line.startswith('#'):
            out_verified.write(mut_line)  # copy header information
        else:
            edit_mutations_by_line(mut_line, gene_locations, frameshifts, out_verified, out_errors, exon_dir)
    mut_inhandle.close()

    out_errors.close()
    out_verified.close()
    sys.stderr.write('Finished updating ' + mut_infile + ' to ' + verified_outfile + '\n' +
                     'See errors in logfile: ' + log_file + '\n')


########################################################################################################
# CREATE MUTATIONAL BURDEN FILES
########################################################################################################

def find_ambiguous_snps(mut_infile, exclude_seqs=None):
    """
    :param mut_infile: full path to a tab-delimited file containing mutation information
    :param exclude_seqs: set of sequence IDs to exclude (pseudogenes?)
    :return: set of SNPs (chromosome, position, alt allele) that are both Silent/Nonsynonymous
    """

    nonsynonymous_snps = set()  # (chrom, position, alt_allele)
    synonymous_snps = set()

    mut_handle = gzip.open(mut_infile) if mut_infile.endswith('gz') else open(mut_infile)
    for mut_line in mut_handle:
        if mut_line.startswith('#'):
            continue

        (sequence_id, _, _, mut_type, genome_build, chrom, chrom_pos,
         ref_allele, alt_allele) = mut_line[:-1].split('\t')[:9]

        if exclude_seqs and sequence_id in exclude_seqs:
            continue

        if mut_type == 'Silent':
            synonymous_snps.add((chrom, chrom_pos, alt_allele))
        else:
            nonsynonymous_snps.add((chrom, chrom_pos, alt_allele))
    mut_handle.close()

    return synonymous_snps.intersection(nonsynonymous_snps)  # ambiguous SNPs


########################################################################################################

def find_expressed_genes(expression_file, limit_expression=True):
    """
    :param expression_file: full path to a tab-delimited file containing genes and the sample IDs where
                           that gene is expressed
    :param limit_expression: boolean indicating whether we want to keep track of samples that express genes or not
    :return: dictionary of gene name -> set(tumor IDs where the gene is expressed)
    """

    # Do we want to filter by expression at all?
    if not limit_expression:
        return None

    expressed_genes = {}

    expr_handle = gzip.open(expression_file) if expression_file.endswith('gz') else open(expression_file)
    for expr_line in expr_handle:
        if expr_line.startswith('#'):
            continue
        gene_name, sample_ids = expr_line[:-1].split('\t')[:2]

        expressed_genes[gene_name.split(',')[0]] = set(sample_ids.split(','))
    expr_handle.close()

    return expressed_genes


########################################################################################################

def get_mutation_counts(verified_mutfile, ambiguous_snps=None, exclude_seqs=None):
    """
    :param verified_mutfile: formatted mutation file of verified mutations
    :param ambiguous_snps: set of silent SNPs that have non-synonymous effects in alternate isoforms
    :param exclude_seqs: set of sequence IDs corresponding to pseudogenes, decayed transcripts, etc. (to ignore)
    :return: dictionary of mutation counts per sample ID AND total number of samples (we assume that each
             mutation file corresponds to one cancer type, though this is not necessarily the case)
    """

    mutations_per_sample = {}

    mutfile_inhandle = gzip.open(verified_mutfile) if verified_mutfile.endswith('.gz') else open(verified_mutfile)
    for mut_line in mutfile_inhandle:
        if mut_line.startswith('#') or len(mut_line[:-1].split('\t')) < 9:
            continue

        sequence_id, _, sample_id, mutation_type, _, chrom_id, chrom_pos, _, alt_allele = mut_line[:-1].split('\t')[:9]

        if sequence_id in exclude_seqs:  # remove pseudogenes ?!
            continue

        if mutation_type == 'Silent' and (chrom_id, chrom_pos, alt_allele) in ambiguous_snps:  # ignore ambiguous silent
            continue

        cancer_types = mut_line[:-1].split('\t')[15] if len(mut_line[:-1].split('\t')) > 15 else ''
        if ('OV.' in verified_mutfile or 'OV' in cancer_types.split(',')) and len(
            sample_id) < 16:  # e.g. TCGA-09-2053-01C
            sample_id = sample_id + 'A'

        if sample_id not in mutations_per_sample:
            mutations_per_sample[sample_id] = {'Nonsynonymous': set()}  # Aggregate count of all nonsynonymous mutations

        if mutation_type not in mutations_per_sample[sample_id]:  # e.g., Missense, Silent, Nonstop, Nonsense
            mutations_per_sample[sample_id][mutation_type] = set()

        mutations_per_sample[sample_id][mutation_type].add((chrom_id, chrom_pos))

        if mutation_type != 'Silent':
            mutations_per_sample[sample_id]['Nonsynonymous'].add((chrom_id, chrom_pos))

    mutfile_inhandle.close()

    return {sample_id: {muttype: len(v.get(muttype, set())) for muttype in v.keys()} for sample_id, v in
            mutations_per_sample.items()}


########################################################################################################

def find_samples_by_cancer(mut_infile, cancer_label):
    """
    :param mut_infile: full path to a formatted mutation file
    :param cancer_label: cancer_type
    :return: sets of sample_ids by cancer_type (i.e., "cancer_label" in the case of non-aggregate files)
    """

    sampleids_by_cancer = {}
    mut_handle = gzip.open(mut_infile) if mut_infile.endswith('.gz') else open(mut_infile)
    for mut_line in mut_handle:
        if mut_line.startswith('#') or len(mut_line[:-1].split('\t')) < 15:
            continue

        sample_id = mut_line[:-1].split('\t')[2]

        cancer_types = mut_line[:-1].split('\t')[15] if len(mut_line[:-1].split('\t')) > 15 else ''
        if ('OV.' in mut_infile or 'OV' in cancer_types.split(',')) and len(sample_id) < 16:  # e.g. TCGA-09-2053-01C
            sample_id = sample_id + 'A'

        if len(cancer_types) > 0:
            cancer_type = [current_cancer for current_cancer in cancer_types.split(',')
                           if current_cancer not in ['GBMLGG', 'COADREAD', 'KIPAN', 'STES', 'Aggregate',
                                                     'AggregateTrunc']][0]
        else:
            cancer_type = cancer_label

        if cancer_type not in sampleids_by_cancer:
            sampleids_by_cancer[cancer_type] = set()
        sampleids_by_cancer[cancer_type].add(sample_id)
    mut_handle.close()

    return sampleids_by_cancer


########################################################################################################

def find_hypermutators(mut_infile, cancer_label, mutations_per_sample, samples_per_cancer, hypermutator_cutoff=2000):
    """
    :param mut_infile: full path to a processed mutation file
    :param cancer_label: cancer_type
    :param mutations_per_sample: dictionary of tumor_sample -> total number of nonsynonymous mutations in that sample
    :param samples_per_cancer: dictionary of cancer_type -> total number of tumor samples with that cancer type
    :param hypermutator_cutoff: minimum number of nonsynonymous mutations per tumor sample to be considered hypermut.
    :return: dictionary of cancer_type -> set(sample_ids) that correspond to "hypermutators," where we define a
             hypermutator as any sample with hypermutator_cutoff nonsynonymous mutations OR has a mutation count
             that is 3 standard deviations above the average mutation count for that cancer type OR has a mutation
             in PolE or PolD1 DNA-repair genes
    """

    hypermutators_per_cancer = {}

    # First, keep track of those samples that have mutations in DNA-repair genes
    mut_handle = gzip.open(mut_infile) if mut_infile.endswith('.gz') else open(mut_infile)
    for mut_line in mut_handle:
        if mut_line.startswith('#') or len(mut_line[:-1].split('\t')) < 15:
            continue

        if 'ENSG00000177084' not in mut_line and 'ENSG00000062822' not in mut_line:  # PolE, PolD1
            continue

        sample_id, muttype = mut_line[:-1].split('\t')[2:4]

        if muttype == 'Silent':
            continue

        cancer_types = mut_line[:-1].split('\t')[15] if len(mut_line[:-1].split('\t')) > 15 else ''
        if ('OV.' in mut_infile or 'OV' in cancer_types.split(',')) and len(sample_id) < 16:  # e.g. TCGA-09-2053-01C
            sample_id = sample_id + 'A'

        if len(cancer_types) > 0:
            cancer_type = [current_cancer for current_cancer in cancer_types.split(',')
                           if current_cancer not in ['GBMLGG', 'COADREAD', 'KIPAN', 'STES', 'PANGI', 'Aggregate',
                                                     'AggregateTrunc']][0]
        else:
            cancer_type = cancer_label

        if cancer_type not in hypermutators_per_cancer:
            hypermutators_per_cancer[cancer_type] = set()
        hypermutators_per_cancer[cancer_type].add(sample_id)

    mut_handle.close()

    # Then, keep track of highly-mutating cancers
    for cancer_type, sampleset in samples_per_cancer.items():
        mutation_distributions = [v.get('Nonsynonymous', 0) for sample_id, v in mutations_per_sample.items() if
                                  sample_id in sampleset]
        mutation_average = arithmetic_mean(mutation_distributions)
        mutation_stddev = standard_deviation(mutation_distributions)
        hypermut_stddev_cutoff = mutation_average + 3 * mutation_stddev

        hypermutators_per_cancer[cancer_type] = set([current_sample for current_sample in sampleset
                                                     if mutations_per_sample[current_sample].get('Nonsynonymous', 0) >=
                                                     min(hypermut_stddev_cutoff, hypermutator_cutoff)])

    return hypermutators_per_cancer


########################################################################################################

def create_mutational_burden_input(mut_infile, mutburden_outfile, expressionfile, cancer_label,
                                   mutation_values=None, exclude_samples=None, exclude_seqs=None,
                                   ambiguous_silent=False):
    """
    :param mut_infile: full path to a reformatted mutation file of verified mutations across all isoforms
    :param mutburden_outfile: outfile to write tab-delimited "mutational burden" output to
    :param expressionfile: full path to a tab-delimited file listing Ensembl gene name and the set of tumor
                           samples, across cancer types, where that gene was expressed at TPM > 0.1
    :param cancer_label: cancer_type
    :param mutation_values: dictionary of mutation -> value (clonal fraction?)
    :param exclude_samples: set of tumor samples to be excluded (hypermutators?)
    :param exclude_seqs: set of genes to be excluded (pseudogenes? decayed?)
    :return: calculate the mutational burden values based on total number of mutations per sample and total
             samples in the mutation file. NOTE that "mutational burden" could also correspond to some outside
             function (i.e., SIFT, PolyPhen2, natural allele frequency, etc.)
    """

    print_errors = False  # print out errors for debugging?
    limit_expression = expression_file is not None  # restrict to mutations in expressed genes?

    # expressed_genes -> {geneID:set(tumorIDs)}
    expressed_genes = {}
    if limit_expression:
        expressed_genes = find_expressed_genes(expressionfile, limit_expression)

    # set((chromosome, chromosome position, alternate allele))
    ambiguous_snps = set()
    if not ambiguous_silent:
        ambiguous_snps = find_ambiguous_snps(mut_infile, exclude_seqs)

    # sample_id -> { mutation_type: count }
    mutations_per_sample = get_mutation_counts(mut_infile, ambiguous_snps, exclude_seqs)

    # sample_grouping (e.g., cancer type) -> set(sampleIDs)
    samples_per_cancer = find_samples_by_cancer(mut_infile, cancer_label)

    # We define "hypermutators" as samples with >2000 mutations OR >mean+3stddev mutations
    # sample_grouping -> mutation_type -> set(sample IDs that are hypermutators)
    hypermutators_per_cancer = find_hypermutators(mut_infile, cancer_label, mutations_per_sample, samples_per_cancer)

    if mutburden_outfile.endswith('.gz'):
        mutburden_outhandle = gzip.open(mutburden_outfile, 'w')
    else:
        mutburden_outhandle = open(mutburden_outfile, 'w')
    mutburden_outhandle.write('# Mutational burden information for ' + mut_infile + '\n')

    allmutations = {}  # sequence_id -> list of hits
    total_mutations = {'Nonsynonymous': set(), 'Synonymous': set()}  # keep track of total mutations in this cancer type

    mutfile_inhandle = gzip.open(mut_infile) if mut_infile.endswith('.gz') else open(mut_infile)
    for mut_line in mutfile_inhandle:
        if mut_line.startswith('#'):
            if '\t' not in mut_line:
                mutburden_outhandle.write(mut_line)
            continue
        if len(mut_line[:-1].split('\t')) < 15:
            if print_errors:
                sys.stderr.write('Improperly formatted line:\n' + mut_line)
            continue

        (sequence_id, other_ids, sample_id, muttype, _, chrom, chrompos, _, altallele,
         _, _, _, prot_pos, ref_aa, alt_aa) = mut_line[:-1].split('\t')[:15]

        # exclude unknown mutations
        if muttype not in ['Missense', 'Silent', 'Nonsense', 'Nonstop']:
            if print_errors:
                sys.stderr.write('Unknown mutation type: ' + muttype + '\n')
            continue

        # exclude ambiguous silent mutations
        if muttype == 'Silent' and (chrom, chrompos, altallele) in ambiguous_snps:
            if print_errors:
                sys.stderr.write('Silent mutation at ' + str(chrom) + ':' + str(chrompos) + ':' +
                                 str(altallele) + ' is ambiguous\n')
            continue

        # update the tumor sample names for OV cancer (in new data)
        cancer_types = mut_line[:-1].split('\t')[15] if len(mut_line[:-1].split('\t')) > 15 else ''
        if ('OV.' in mut_infile or 'OV' in cancer_types.split(',')) and len(sample_id) < 16:  # e.g. TCGA-09-2053-01C
            sample_id = sample_id + 'A'

        # exclude genes that are unexpressed in particular samples?
        ensembl_id = other_ids[other_ids.find('gene_id:') + len('gene_id:'):].split(',')[0]
        if expressed_genes and ensembl_id not in expressed_genes:
            if print_errors:
                sys.stderr.write('Could not find ' + ensembl_id + ' in expressed_genes.\n')
            continue

        if expressed_genes and ensembl_id in expressed_genes and sample_id not in expressed_genes[ensembl_id]:
            if print_errors:
                sys.stderr.write('Could not find ' + sample_id + ' in expressed_genes[' + ensembl_id + ']\n')
            continue

        # exclude hypermutators
        if exclude_samples and sample_id in exclude_samples:
            if print_errors:
                sys.stderr.write(sample_id + ' was found in exclude_samples\n')
            continue

        # exclude pseudogenes / decayed / truncated genes
        if exclude_seqs and sequence_id in exclude_seqs:
            if print_errors:
                sys.stderr.write(sequence_id + ' was found in exclude_seqs\n')
            continue

        # If we are given an aggregate file, exclude_samples and exclude_seqs become cancer-specific
        #  dictionaries (e.g., {'ACC'->set(), 'BLCA'->set()})
        # aggregate file case
        if len(cancer_types) > 0:
            skip = False
            for file_type in [cancer_label] + cancer_types.split(','):
                if (exclude_samples and file_type in exclude_samples and sample_id in exclude_samples[file_type]) or \
                   (exclude_seqs and file_type in exclude_seqs and sequence_id in exclude_seqs[file_type]):
                    skip = True
            if skip:
                continue

            # reset the file label to non-aggregate sets...
            cancer_label = [current_cancer for current_cancer in cancer_types.split(',')
                            if current_cancer not in ['GBMLGG', 'COADREAD', 'KIPAN', 'STES', 'PANGI', 'Aggregate',
                                                      'AggregateTrunc']][0]

        # Add the mutations for this sequence (and keep track of total mutations) !!!
        subclonal_frac = 1 if not mutation_values else mutation_values.get((sample_id, chrom, chrompos), 1)
        if sequence_id not in allmutations:
            allmutations[sequence_id] = []
        total_mutations['Synonymous' if muttype == 'Silent' else 'Nonsynonymous'].add(
            (sample_id, chrom, chrompos, subclonal_frac))
        allmutations[sequence_id].append((int(prot_pos) - 1,  # 0-index protein position
                                          muttype,  # Missense or Silent?
                                          subclonal_frac,  # subclonal fraction
                                          mutations_per_sample[sample_id][muttype],  # total mutations in sample
                                          sample_id + ',' + str(sample_id in hypermutators_per_cancer[cancer_label]),
                                          len(samples_per_cancer[cancer_label]),  # number tumor samples per cancer
                                          ref_aa,  # reference amino acid
                                          alt_aa,  # alternate amino acid
                                          cancer_label,  # cancer type
                                          sample_id))  # sample ID (ensure otherwise identical entries are not excluded)
    mutfile_inhandle.close()

    for muttype in total_mutations.keys():
        mutburden_outhandle.write('# Total ' + muttype + ' mutations = ' + str(len(total_mutations[muttype])) +
                                  ' | ' + str(len([sample_id for sample_id, _, _, _ in total_mutations[muttype]
                                                   if
                                                   sample_id not in hypermutators_per_cancer[cancer_label]])) + '\n' +
                                  '# Total ' + muttype + ' mutation values = ' +
                                  str(sum([mutburden for _, _, _, mutburden in total_mutations[muttype]])) +
                                  ' | ' + str(sum([mutburden for sample_id, _, _, mutburden in total_mutations[muttype]
                                                   if
                                                   sample_id not in hypermutators_per_cancer[cancer_label]])) + '\n' +
                                  '# Total squared ' + muttype + ' mutation values = ' +
                                  str(sum([mutburden ** 2 for _, _, _, mutburden in total_mutations[muttype]])) +
                                  ' | ' + str(
            sum([mutburden ** 2 for sample_id, _, _, mutburden in total_mutations[muttype]
                 if sample_id not in hypermutators_per_cancer[cancer_label]])) + '\n')
    mutburden_outhandle.write('\t'.join(['#SequenceID', 'Sequence-0-Index', 'MutationType', 'MutationValue',
                                         'MutationsPerSample', 'SampleID,Hypermutator?', 'SamplesPerCancer',
                                         'ReferenceAminoAcid', 'AlternateAminoAcid', 'CancerType']) + '\n')
    for sequence_id in sorted(allmutations.keys()):
        for m in sorted(allmutations[sequence_id]):
            # Sample ID was added just to ensure otherwise identical entries would not be overwritten.
            mutburden_outhandle.write(sequence_id + '\t' + '\t'.join(map(str, m[:-1])) + '\n')
    mutburden_outhandle.close()
    sys.stderr.write('Mutational burden file in ' + mutburden_outfile + '\n')


########################################################################################################

def get_mutation_clonal_fraction(maf_file):
    """
    :param maf_file: full path to a MuSE (or other GDC) maf file
    :return: dictionary of (tumor sample ID, chromosome, location) -> subclonal fraction
    """

    mutation_to_clonalfrac = {}

    if not os.path.isfile(maf_file):
        sys.stderr.write('Could not open ' + maf_file + '\n')
        sys.exit(1)

    maf_handle = gzip.open(maf_file) if maf_file.endswith('gz') else open(maf_file)
    header = None
    for maf_line in maf_handle:
        if maf_line.startswith('#'):
            continue
        header = maf_line.lower()[:-1].split('\t')
        break

    for maf_line in maf_handle:
        m = maf_line[:-1].split('\t')
        sample_id = '-'.join(m[header.index('tumor_sample_barcode')].split('-')[:4])
        if 'OV.' in maf_file and len(sample_id) < 16:  # e.g. TCGA-09-2053-01C
            sample_id = sample_id + 'A'

        chromosome = m[header.index('chromosome')]
        if chromosome.startswith('chr'):
            chromosome = chromosome[3:]
        chrom_position = m[header.index('start_position')].split('-')[0]
        mutational_burden = float(m[header.index('t_alt_count')]) / float(m[header.index('t_depth')])

        mutation_to_clonalfrac[(sample_id, chromosome, chrom_position)] = mutational_burden
    maf_handle.close()

    return mutation_to_clonalfrac


########################################################################################################
# MAIN
########################################################################################################

if __name__ == "__main__":

    # Change any of the following "if False:" prompts to "if True:" in order to run the relevant code

    parser = argparse.ArgumentParser(description='Properly format cancer mutation files to run.')

    parser.add_argument('--cancer_type', type=str,
                        help='Subset of mutations to run on from a particular cancer subtype (e.g., Aggregate)',
                        default='Aggregate',
                        choices=['ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC', 'KICH',
                                 'KIRC', 'KIRP', 'LAML', 'LCML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV', 'PAAD',
                                 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA', 'THYM', 'UCEC', 'UCS',
                                 'UVM', 'COADREAD', 'GBMLGG', 'KIPAN', 'PANGI', 'STES', 'AggregateTrunc', 'Aggregate'])
    parser.add_argument('--bias_type', type=str, default='both', choices={'both', 'missense', 'gc'},
                        help='Type of mutational bias to include when creating per-protein-position mut. bias tracks.')

    parser.add_argument('--mutation_path', type=str, default=data_path + 'gdc/somatic_mutations/',
                        help='Full path to a directory where .maf files are stored (and where output files will also' +
                             ' be stored); subdirectories correspond to cancer types')
    parser.add_argument('--exon_path', type=str, default=data_path + 'ensembl/Homo_sapiens.GRCh38/exons/',
                        help='Full path to a directory where subdirectories contain genes and per-protein exon info.')

    parser.add_argument('--expression_file', type=str,
                        default=data_path + 'gdc/expression/TCGA_GRCh38_expressed-genes_TPM.tsv.gz',
                        help='Full path to a tab-delimited list of genes and which tumor samples they are expressed ' +
                             'in at >0.1 TPM')
    parser.add_argument('--seq_file', type=str,
                        default=data_path + 'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.pep.all.withgenelocs.' +
                                'verified.fa.gz',
                        help='Full path to a fasta-formatted sequence file where headers contain exon location ' +
                             'information as well as transcript status (decay, etc.)')
    parser.add_argument('--cdna_file', type=str,
                        default=data_path + 'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.cdna.all.withgenelocs.' +
                                'verified.fa.gz',
                        help='Full path to a fasta-formatted file containing cDNA sequences for each protein isoform.')
    parser.add_argument('--mutbias_file', type=str,
                        default=data_path+'gdc/somatic_mutations/TCGA_GRCh38_mutational-biases.tsv.gz',
                        help='Full path to a tab-delimited list of all the C/G mutation biases per cancer.')
    parser.add_argument('--manifest', type=str,
                        default=data_path + 'gdc/somatic_mutations/gdc_manifest.muse-masked_somatic_mutation.tsv',
                        help='Full path to a GDC manifest file containing directory and file name information.')

    parser.add_argument('--rename_files', dest='rename_files', action='store_true', default=False,
                        help='Rename files downloaded from GDC to something more intuitive to work with.')
    parser.add_argument('--aggregate_files', dest='aggregate_files', action='store_true', default=False,
                        help='Create aggregate .maf files in the appropriate directories.')
    parser.add_argument('--verify_mutations', dest='verify_mutations', action='store_true', default=False,
                        help='Step 1: Verify that mutations are valid (reference matches what is expected)')
    parser.add_argument('--create_mut_burden', dest='create_mut_burden', action='store_true', default=False,
                        help='Step 2: Create mutational burden files from verified mutation files')
    parser.add_argument('--mutational_bias', dest='mutational_bias', action='store_true', default=False,
                        help='Create a file specifying missense and GC biases across mutations per cancer type.')
    parser.add_argument('--create_mutbias_tracks', dest='create_mutbias_tracks', action='store_true', default=False,
                        help='Create a file of per-protein-position mutational biases for a specific cancer type.')

    parser.add_argument('--include_unexpressed', dest='include_unexpressed', action='store_true', default=False,
                        help='Include mutations from genes that were not expressed at TPM > 0.1 in the sample?')
    parser.add_argument('--include_decay', dest='include_decay', action='store_true', default=False,
                        help='Include mutations from those "proteins" that resulted in transcript decay?')
    parser.add_argument('--include_frameshift', dest='include_frameshift', action='store_true', default=False,
                        help='Include mutations that occur after a frameshift or nonsense mutation?')
    parser.add_argument('--include_ambiguous', dest='include_ambiguous', action='store_true', default=False,
                        help='Include silent mutations with nonsynonymous effects in other protein isoforms?')

    args = parser.parse_args()

    # ---------------------------------------------------------------------------------------------------------
    if not os.path.isdir(args.mutation_path):
        sys.stderr.write('Could not find mutation path in ' + args.mutation_path + '\n' +
                         'Please see https://github.com/Singh-Lab/pertinint-internal/wiki/Cancer-mutation-data for ' +
                         'instructions on downloading the appropriate mutation files to this directory.\n')
        sys.exit(1)

    # ---------------------------------------------------------------------------------------------------------
    if args.rename_files:
        """
        Rename the files automatically downloaded from NCI's Genomic Data Commons to intuitive names
        """

        if not os.path.isfile(args.manifest):
            sys.stderr.write('Could not open manifest file: ' + args.manifest + '\n' +
                             'Please run with: python ' + sys.argv[0] + ' --manifest <manifest_file> --rename_files\n')
            sys.exit(1)

        cancer_types = {
            'ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'DLBC', 'ESCA', 'GBM', 'HNSC',
            'KICH', 'KIRC', 'KIRP', 'LAML', 'LGG', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV',
            'PAAD', 'PCPG', 'PRAD', 'READ', 'SARC', 'SKCM', 'STAD', 'TGCT', 'THCA',
            'THYM', 'UCEC', 'UCS', 'UVM', 'COADREAD', 'GBMLGG', 'KIPAN', 'PANGI', 'STES'}

        with open(args.manifest) as manifest_handle:
            manifest_handle.next()  # skip header
            for manifest_line in manifest_handle:
                directory, filename = manifest_line[:-1].split('\t')[:2]
                for c_type in cancer_types:
                    if '.' + c_type + '.' in filename:
                        current_cancer_type = c_type
                        break
                else:
                    continue

                # create a new directory and move the appropriate file there
                if not os.path.isdir(args.mutation_path + current_cancer_type):
                    call(['mkdir', args.mutation_path + current_cancer_type])
                if os.path.isdir(args.mutation_path + directory):
                    call(['mv', args.mutation_path + directory + '/annotations.txt',
                          args.mutation_path + current_cancer_type + '/'])
                    call(['mv', args.mutation_path + directory + '/' + filename,
                          args.mutation_path + current_cancer_type + '/'])
                    call(['rm', '-rf', args.mutation_path + directory])

    # ---------------------------------------------------------------------------------------------------------
    if args.aggregate_files:
        """
        Create aggregate .maf files for specific cancer types ('COADREAD', 'GBMLGG', 'KIPAN', 'PANGI', 'STES',
        'Aggregate', 'AggregateTrunc')
        """

        agg_groups = aggregate_groups()
        for agg_cancer_type, cancer_subtypes in agg_groups.items():

            # create output aggregated .maf file
            if not os.path.isdir(args.mutation_path + agg_cancer_type):
                call(['mkdir', args.mutation_path + agg_cancer_type])

            output_maf = args.mutation_path+agg_cancer_type+'/TCGA.'+agg_cancer_type+'.muse.aggregated.somatic.maf.gz'

            # get list of input .maf files to be included
            input_mafs = set()
            for cancer_type in [ctype for ctype in os.listdir(args.mutation_path) if ctype in cancer_subtypes]:
                mutation_file = get_maf_file(args.mutation_path, cancer_type, False)
                if mutation_file:
                    input_mafs.add(mutation_file)

            # create the aggregate maf file
            if create_aggregate_maf_files(input_mafs, output_maf):
                sys.stderr.write('Wrote to ' + output_maf + '\n')

    # ---------------------------------------------------------------------------------------------------------
    if args.verify_mutations:
        """
        Process an input .maf file to create two new intermediate files: 
        (1) .exon_mutation_list.txt = subset of mutations landing in exon regions in genes that have been verified
        (2) .verified_mutation_list.txt = predicted functional consequences of the exonic mutations based on 
            mapping to verified cDNA sequences and translating to protein
        """

        # (1) check that the required arguments are valid:
        if not os.path.isdir(args.exon_path):
            sys.stderr.write('Could not find exon path in ' + args.exon_path + '\n' +
                             'Please see https://github.com/Singh-Lab/pertinint-internal/wiki/' +
                             'Preprocess-protein-sequences for instructions on\ncreating this directory and ' +
                             'generating the appropriate input files.\n')
            sys.exit(1)
        if not os.path.isfile(args.seq_file):
            sys.stderr.write('Could not find FASTA sequences in ' + args.seq_file + ', ' +
                             'set variable with --seq_file flag.\n')
            sys.exit(1)

        maf_file = get_maf_file(args.mutation_path, args.cancer_type, True)  # maf file that we are interested in
        sys.stderr.write('Processing .maf file: ' + maf_file + '\n')

        # (2) extract exonic mutations that should later be verified with respect to all protein isoforms:
        exon_file = args.mutation_path+args.cancer_type+'/TCGA.'+args.cancer_type+'.muse.processed-exonic.somatic.tsv'
        extract_mutations_from_maf(maf_file, exon_file)

        # (3) verify exonic mutations with respect to all protein isoforms:
        if os.path.isfile(exon_file):
            mut_file = args.mutation_path + args.cancer_type + '/TCGA.' + args.cancer_type + \
                       '.muse.processed-verified.somatic.tsv'
            remap_mutations_to_all_proteins(exon_file,  # full path to reformatted mutations occurring in verified exons
                                            mut_file,  # full path to an output file to write verified mutations to
                                            mut_file.replace('.tsv', '-ERRORS.txt'),  # output file to write errors to
                                            args.seq_file,  # fasta file where headers include relevant information
                                            args.exon_path,  # path to directory with gene/protein/exon information
                                            args.include_frameshift)  # include mutations after frameshifts?

            sys.stderr.write('Verified mutations found in ' + mut_file + '\n')
        else:
            sys.stderr.write('Failed to find any exonic mutations in ' + maf_file + '/failed to write to ' + exon_file +
                             '...Exiting.\n')
            sys.exit(1)

    # ---------------------------------------------------------------------------------------------------------
    if args.create_mut_burden:
        """
        Create new "mutational burden" files from the verified mutation data that identifies hypermutating samples,
        restricts to mutations in proteins that are expressed at TPM > 0.1 and whose transcripts are not "decayed"
        """

        # (1) check if our verified mutation information and .maf file are both available:
        if not os.path.isfile(args.seq_file):
            sys.stderr.write(
                'Could not find FASTA sequences in ' + args.seq_file + ', set variable with --seq_file flag.\n')
            sys.exit(1)

        mut_file = args.mutation_path + args.cancer_type + '/TCGA.' + args.cancer_type + '.muse.processed-verified.' + \
                   'somatic.tsv'  # verified mutation file
        if not os.path.isfile(mut_file):
            sys.stderr.write('Failed to find verified mutations in ' + mut_file + '\n' +
                             'Run: python ' + sys.argv[
                                 0] + ' --verify_mutations --cancer_type ' + args.cancer_type + '\n')
            sys.exit(1)

        maf_file = get_maf_file(args.mutation_path, args.cancer_type, True)  # input maf file

        # (2) Exclude proteins that are pseudogenes or marked with "decay":
        exclude_seqids = set()
        if not args.include_decay:
            sys.stderr.write('Finding pseudogenes and decayed transcripts to exclude...\n')
            ignore_flags = ['transcript_biotype:' + fault for fault in ['non_stop_decay', 'nonsense_mediated_decay',
                                                                        'polymorphic_pseudogene']]
            fasta_handle = gzip.open(args.seq_file) if args.seq_file.endswith('gz') else open(args.seq_file)
            for seq_line in fasta_handle:
                if seq_line.startswith('>') and True in [flag in seq_line for flag in ignore_flags]:
                    exclude_seqids.add(seq_line[1:-1].split()[0])
            fasta_handle.close()
            sys.stderr.write('Done!\n')

        # (3) Specify file with expression data if we are excluding lowly/unexpressed genes:
        expression_file = None
        if not args.include_unexpressed:
            if not os.path.isfile(args.expression_file):
                sys.stderr.write(
                    'Could not find tumor sample expression information in ' + args.expression_file + '\n' +
                    'See https://github.com/Singh-Lab/pertinint-internal/wiki/Cancer-mutation-data ' +
                    'for instructions on generating this file.\n')
                sys.exit(1)
            expression_file = args.expression_file

        # (4) exclude particular tumor samples if desired:
        exclude_tumors = set()

        # (5) Extract "mutational burden" (clonal fraction) values from the specified .maf file:
        mutation_values = get_mutation_clonal_fraction(maf_file)  # dictionary of sample ID -> something.

        # (6) Finally, create the mutational burden file!
        out_mutburden_file = args.mutation_path + args.cancer_type + '/TCGA.' + args.cancer_type + \
                             '.muse.processed-mutational-burden.somatic.tsv'  # mutational burden file!
        create_mutational_burden_input(mut_file,  # verified mutation file
                                       out_mutburden_file,  # output file to write to
                                       expression_file,  # file with lists of expressed genes by tumor sample
                                       args.cancer_type,  # cancer type to run on
                                       mutation_values,  # dictionary of mutation -> clonal fraction
                                       exclude_tumors,  # specific tumor samples to exclude
                                       exclude_seqids,  # specific protein sequences to exclude
                                       args.include_ambiguous)  # include silent mutations with nonsilent effects?

    # ---------------------------------------------------------------------------------------------------------
    if args.mutational_bias:
        """
        Create a mutational bias file to be used when preprocessing tracks
        """

        # get list of input .maf files to be included
        input_mafs = {}
        for cancer_type in [ct for ct in os.listdir(args.mutation_path) if os.path.isdir(args.mutation_path+ct)]:
            mutation_file = get_maf_file(args.mutation_path, cancer_type, False)
            if mutation_file:
                if cancer_type not in input_mafs:
                    input_mafs[cancer_type] = []
                input_mafs[cancer_type].append(mutation_file)

        # create the mutational bias file:
        mutation_bias_file = args.mutbias_file
        if summarize_mutations_by_cancer(input_mafs, mutation_bias_file):
            sys.stderr.write('Wrote to ' + mutation_bias_file + '\n')

    # ---------------------------------------------------------------------------------------------------------
    if args.create_mutbias_tracks:
        """
        Create a (large) file containing per-protein-position likelihoods of harboring a missense mutation
        """

        if not os.path.isfile(args.cdna_file):
            sys.stderr.write(
                'Could not find FASTA cDNA sequences in ' + args.seq_file + ', ' +
                'set variable with --cdna_file flag.\n' +
                'Please run: python verify_sequences.py --create_final_fasta')
            sys.exit(1)

        if args.bias_type in ['both', 'gc'] and not os.path.isfile(args.mutbias_file):
            sys.stderr.write(
                'Could not find list of per-cancer mutational biases in '+args.mutbias_file+'\n' +
                'Please run: python ' + sys.argv[0] + ' --mutational_bias --mutbias_file '+args.mutbias_file+'\n')
            sys.exit(1)

        final_frequency_file = args.mutation_path + args.cancer_type + '/TCGA.' + args.cancer_type + \
                               '.muse.mutational_biases.tsv'
        if get_mutation_likelihood(args.cancer_type,
                                   args.cdna_file,
                                   args.mutbias_file,
                                   final_frequency_file,
                                   args.bias_type):
            sys.stderr.write('Wrote to ' + final_frequency_file + '\n')
