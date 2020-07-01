#!/usr/bin/python

"""
Create a BLAST-P database from a FASTA file containing all PDB chain sequences that ConCavity
ran on

e.g.,
>101m_A
MADFHGJAJTHAR..

>101m_B
ADFGKAHETG

Then, run blastp on a NEW fasta file (all Ensembl human protein sequences) on that newly created
database to get matches
"""

import os
import sys
import gzip
from subprocess import call
from config import data_path
from process_weight_vectors import get_sequence_lengths, get_original_background_mutation, get_background_variation, update_domain_weights, empirically_determine_minmuts


########################################################################################################
# CONSTANTS
########################################################################################################

BUILD = 'GRCh38'
CANCER_TYPE = 'Aggregate'


########################################################################################################

def run_blastp_search_against_concavity(dir_to_blastdb, pdb_fasta, query_fasta, blast_results):
    """
    :param dir_to_blastdb: full path to a directory where you want the formatted BLAST-P database
                           files to be stored (your PDB fasta file should be the only file in here)
    :param pdb_fasta: full path to a fasta-formatted file containing the protein sequences of all the PDB
                      chains that ConCavity ran on (should be in the dir_to_blastdb directory)
    :param query_fasta: full path to an UNZIPPED fasta file containing all human protein sequences that
                        you want to query against your newly created blast database
    :param blast_results: full path to a file where you want to store the BLAST output
    :return: None, but print success messages at each step in the process.
    """

    # first, confirm that your input pdb_fasta file exists
    if not os.path.isfile(pdb_fasta):
        sys.stderr.write('Could not find PDB fasta file: ' + pdb_fasta + '\n')
        sys.exit(1)

    # then, confirm that your blast directory exists (create otherwise)
    if dir_to_blastdb.endswith('/'):
        dir_to_blastdb = dir_to_blastdb[:-1]

    if not os.path.isdir(dir_to_blastdb):
        sys.stderr.write('Attempting to create ' + dir_to_blastdb + '...\n')
        call(['mkdir', dir_to_blastdb])
    if not os.path.isdir(dir_to_blastdb):
        sys.stderr.write('Could not create ' + dir_to_blastdb + '\n')
        sys.exit(1)

    # move your fasta file into your BLAST directory if it's not already in there:
    if pdb_fasta != dir_to_blastdb + '/' + pdb_fasta.split('/')[-1]:
        sys.stderr.write('Copying ' + pdb_fasta + ' to ' + dir_to_blastdb + '/...\n')
        call(['cp', pdb_fasta, dir_to_blastdb + '/'])
        pdb_fasta = dir_to_blastdb + '/' + pdb_fasta.split('/')[-1]

    # unzip pdb fasta file if need be...
    if pdb_fasta.endswith('.gz'):
        sys.stderr.write('Attempting to unzip ' + pdb_fasta + '...\n')
        call(['gzip', '-d', pdb_fasta])
    if not os.path.isfile(pdb_fasta.replace('.gz', '')):
        sys.stderr.write('Could not unzip ' + pdb_fasta + '\n')
        sys.exit(1)

    # create the blastp database!
    system_call = ['makeblastdb', '-in', pdb_fasta, '-out', dir_to_blastdb]
    sys.stderr.write('Creating blastp database: ' + ' '.join(system_call) + '...\n')
    call(system_call)
    sys.stderr.write('Success!\n')

    # confirm existence of query fasta file (and unzip if necessary)
    if not os.path.isfile(query_fasta):
        sys.stderr.write('Could not find query fasta file: ' + query_fasta)
        sys.exit(1)
    if query_fasta.endswith('.gz'):
        sys.stderr.write('Attempting to unzip ' + query_fasta + '...\n')
        call(['gzip', '-d', query_fasta])
    if not os.path.isfile(query_fasta.replace('.gz', '')):
        sys.stderr.write('Could not unzip ' + query_fasta + '\n')
        sys.exit(1)

    # finally, run blastp !
    sys.stderr.write('Running blastp, storing results in ' + blast_results + '\n')
    system_call = ['blastp', '-query', query_fasta.replace('.gz', ''),
                   '-db', dir_to_blastdb + '/' + query_fasta.split('/')[-1].split('.')[0],
                   '-evalue', '1E-6',
                   '-outfmt', '6',
                   '-num_threads', '4',
                   '-out', blast_results]
    call(system_call)
    sys.stderr.write('Success!\n')


########################################################################################################

def process_concavity_results(concavity_scores=data_path + 'concavity/concavity_final_output_032118.txt',
                              structural_coverage=data_path + 'concavity/concavity_bitstring_output_032118.txt'):
    """
    :return:
    """

    sys.stderr.write('Loading sequence/mutation rate information...\n')
    sequence_lengths, sequence_locs = get_sequence_lengths()
    all_mutation_rates = get_original_background_mutation(data_path + 'tcga/mutational-biases-per-protein_' +
                                                          BUILD + '_' + CANCER_TYPE + '.txt.gz',
                                                          set(sequence_lengths.keys()))
    sys.stderr.write('Done!\n')

    if structural_coverage:
        struct_handle = gzip.open(structural_coverage) if structural_coverage.endswith('gz') else open(
            structural_coverage)

    header = '\n'.join(['# ConCavity (v0.1) was downloaded from http://compbio.cs.princeton.edu/concavity/' +
                        'concavity_distr.tar.gz',
                        '#  and run on 136,162 structures downloaded from the PDB on January 4, 2018.',
                        '# ConCavity was run *without* sequence conservation information',
                        '\t'.join(['# Protein_ID', 'Track_ID (of 1)', 'Track_Name', '0-index_Enrichment_Interval(s)',
                                   '0-index_Non-Zero_Binding_Weights', 'Expectation_yi', 'Variance_yi', 'Covariance',
                                   'Minimum_Mutation_Count', 'Empirical_Runtime'])])

    concav_handle = gzip.open(concavity_scores) if concavity_scores.endswith('gz') else open(concavity_scores)
    for concav_line in concav_handle:
        prot_id, scores = concav_line.split()

        # model multiple modeled binding pockets
        if structural_coverage:  # we must increment the bitstring file handle also so that they are synced up!
            struct_coverage = map(bool, map(int, list(struct_handle.next().split()[1])))

        if prot_id not in sequence_locs:
            continue

        if not structural_coverage:
            struct_coverage = [True] * sequence_lengths[prot_id]

        score_dict = {int(ind.split(':')[0]): ind.split(':')[1] for ind in scores.split(',') if
                      float(ind.split(':')[1]) > 0}
        chrom, gene = sequence_locs[prot_id]
        outfile = data_path + 'weightvectors/GRCh38_nucleotide-bias/Aggregate/' + chrom + '/' + gene + '/' + prot_id + \
                  '.concavity-weightvector-orig.tsv.gz'
        outhandle = gzip.open(outfile, 'w') if outfile.endswith('gz') else open(outfile, 'w')
        outhandle.write(header + '\n')

        intervals = []
        start, end = None, None
        for i, cov in enumerate(struct_coverage):
            if cov:
                if not start:
                    start = i
                end = i  # update end each time
            else:  # if we reach False
                if start:  # add the previous interval, if necessary
                    intervals.append((start, end))
                start = None  # reset until the next True
                end = None
        if start and (len(intervals) < 1 or intervals[-1][0] != start):
            intervals.append((start, end))

        mutation_rates = {index: variation for index, variation in
                          enumerate(get_background_variation(all_mutation_rates, prot_id, False))}

        for i, (start, end) in enumerate(intervals):
            outhandle.write(
                '\t'.join([prot_id, str(i + 1), 'ConCavity_' + str(start) + '-' + str(end), str(start) + '-' + str(end),
                           ','.join([str(ind) + ':' + score_dict[ind] for ind in range(start, end + 1)
                                     if ind in score_dict])]) + '\n')
        outhandle.close()
        sys.stderr.write('Wrote to ' + outfile + '\n')

        update_domain_weights(outfile, outfile.replace('-orig', '-exp'), mutation_rates)
        call(['rm', outfile])
        empirically_determine_minmuts(outfile.replace('-orig', '-exp'), outfile.replace('-orig', ''), mutation_rates)
        call(['rm', outfile.replace('-orig', '-exp')])
    concav_handle.close()
    if structural_coverage:
        struct_handle.close()


########################################################################################################
# MAIN
########################################################################################################

if __name__ == "__main__":
    # full path to where you want your BLAST DATABASE to be stored
    blast_db = '/home/snadimpa/datadb/concavity/blastdb'

    # full path to the FASTA file you created from the ConCavity results
    struc_fasta = '/home/snadimpa/datadb/concavity/blastdb/concavity_results.fa'

    # full path to the Ensembl protein sequences fasta file
    prot_fasta = '/home/snadimpa/datadb/ensembl/Homo_sapiens.GRCh38.pep.all.verified.fa'

    # full path to where you want your BLAST-P results to be stored
    blast_output = '/home/snadimpa/datadb/concavity/pdb_to_ensembl_matches.tsv'

    # run blast p !
    run_blastp_search_against_concavity(blast_db, struc_fasta, prot_fasta, blast_output)

    # some code that really belongs in process_weight_vectors.py:
    # process_concavity_results()
