#!/usr/bin/python

"""
Given formatted by-HMM HMMER hits, concatenate all results to remove duplications and order results, 
then also create a file listing ALL complete domains that passed their gathering thresholds.

Contact snadimpa@princeton.edu with questions.
"""

import os
import sys
import gzip
import math
import argparse

####################################################################################################
# CONSTANTS
####################################################################################################

# Release dates and number of entries can always be found for all Pfam releases at
# ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt
HMMINFO = {'27': {'date': 'March 2013', 'entries': '14,831'},
           '28': {'date': 'May 2015', 'entries': '16,230'},
           '29': {'date': 'December 2015', 'entries': '16,295'},
           '30': {'date': 'June 2016', 'entries': '16,306'},
           '31': {'date': 'March 2017', 'entries': '16,712'},
           '32': {'date': 'October 2018', 'entries': '17,929'}}

# Full path to the directory containing this script:
REPO_DIR = os.path.dirname(os.path.abspath(__file__)) + '/'


####################################################################################################

def create_allhmmresbyprot(fasta_infile,
                           pfam_version='32',
                           results_directory=REPO_DIR + 'domains/processed-v32/',
                           outfile=REPO_DIR + 'domains/allhmmresbyprot-v32.tsv.gz'):
    """
    :param fasta_infile: full path to an input FASTA file that HMMER was run on (for header information)
    :param pfam_version: Which Pfam version to use? 32 is our standard
    :param results_directory: location of all .hmmres.gz files generated from the processhmmer.py script
    :param outfile: final file of all non-redundant domain hits
    :return: concatenate all HMMER results into a single file, removing duplicates as best we can
             between the two versions of HMMER
    """
    output_handle = gzip.open(outfile, 'w') if outfile.endswith('gz') else open(outfile, 'w')
    output_handle.write(
        '# HMMER 2.3.2 and HMMER 3.1b2 results on all protein sequences found in ' + fasta_infile + '\n')
    output_handle.write('# Pfam version ' + pfam_version + '.0, released ' + HMMINFO[pfam_version]['date'] +
                        ', containing ' + HMMINFO[pfam_version]['entries'] + ' entries\n')

    output_handle.write(
        '\t'.join(['#TargetID', 'HMM_Name', 'E-value', 'BitScore', 'TargetStart', 'TargetEnd', 'HMM_Seq',
                   'Target_Seq', 'HMM_Pos', 'Description']) + '\n')

    total_domain_hits = 0  # total number of domains processed

    for fname in sorted([a for a in os.listdir(results_directory)
                         if a.startswith('PF') and a.endswith('-v' + pfam_version + '.hmmres.gz')]):

        # (targetID, start, end, pfamHMM, desc) -> (bitscore, evalue, HMM sequence, target sequence, HMM match states)
        results = {}

        results_handle = gzip.open(results_directory + fname)
        for res_line in results_handle:
            if not res_line.startswith('#') and len(res_line.strip()) > 0:
                targ_id, hmm_id, evalue, bitscore, tstart, tend, hmmseq, tseq, hmmpos, desc = res_line[:-1].split('\t')

                # This is a new result; add it:
                if (targ_id, int(tstart), int(tend), hmm_id, desc) not in results:
                    results[(targ_id, int(tstart), int(tend), hmm_id, desc)] = [float(bitscore), evalue, hmmseq, tseq,
                                                                                hmmpos]

                # This is a DUPLICATE result; update if it is better:
                if float(bitscore) > results[(targ_id, int(tstart), int(tend), hmm_id, desc)][0]:
                    results[(targ_id, int(tstart), int(tend), hmm_id, desc)] = [float(bitscore), evalue, hmmseq, tseq,
                                                                                hmmpos]
        results_handle.close()

        # Write out all sorted results:
        for ((targ_id, tstart, tend, hmm_id, desc),
             (bitscore, evalue, hmmseq, tseq, hmmpos)) in sorted(results.items()):
            output_handle.write('\t'.join([targ_id, hmm_id, evalue, str(bitscore), str(tstart), str(tend),
                                           hmmseq, tseq, hmmpos, desc]) + '\n')

        total_domain_hits += len(results.keys())  # keep track of grand total across all domain types
    output_handle.close()

    sys.stderr.write('Concatenated all files from ' + results_directory + ' to ' + outfile + '\n')
    sys.stderr.write(str(total_domain_hits) + ' total domains!\n')


####################################################################################################

def find_domains_from_file(concatenated_file=REPO_DIR + 'domains/processed-v32/allhmmresbyprot-v32.tsv.gz'):
    """
    :param concatenated_file: full path to a formatted HMMER result file
    :return: set of all domains found in that HMMER result file
    """

    dom_handle = gzip.open(concatenated_file) if concatenated_file.endswith('.gz') else open(concatenated_file)
    allhmms = set([a.strip().split('\t')[1] for a in dom_handle if not a.startswith('#')])
    dom_handle.close()

    return allhmms


####################################################################################################

def get_high_information_content(hmmfile):
    """
    :param hmmfile: full path to a Pfam HMM file
    :return: a dictionary of 1-indexed match state -> required amino acid assignment if the
       information content at that match state is >=4 (corresponding to ~95% the same amino acid at
       that position)
    """

    requiredstates = {}

    reach_hmm = False  # boolean indicating whether the match state probabilities have been reached in the HMM file
    aas = []
    with open(hmmfile) as hmm_handle:
        for hmm_line in hmm_handle:
            if hmm_line.startswith('HMM') and not hmm_line.startswith('HMMER'):
                aas = hmm_line.strip().split()[1:]  # amino acids in order
                hmm_handle.next()  # description of transition probabilities
                hmm_handle.next()  # begin state match state emission probabilities (unnecessary)
                hmm_handle.next()  # begin state insertion state emission probabilities (also unnecessary)
                hmm_handle.next()  # begin state transition probabilities
                reach_hmm = True  # Have we reached the actual HMM model yet?
                continue
            elif reach_hmm and len(hmm_line.strip().split()) > len(aas):
                matchstate = hmm_line.strip().split()[0]
                probabilities = map(lambda var: math.exp(-1 * float(var)), hmm_line.strip().split()[1:len(aas) + 1])
                if math.log(20, 2) + sum(j * math.log(j, 2) for j in probabilities) >= 4:  # required amino acid!
                    requiredstates[matchstate] = sorted(zip(probabilities, aas), reverse=True)[0][1]
    return requiredstates


####################################################################################################

def return_passing_hits(current_hits, gathering_thresholds):
    """
    In the case of repetitive domains, multiple domain hits is evidence of biological importance,
    so domains are included even if individually they did not pass the gathering threshold:
    :param current_hits: hmmID -> start_position -> (tuple(hmm_match_states), bitscore)
    :param gathering_thresholds: hmmID -> (sequence gathering threshold, domain gathering threshold)
    :return: the set of HMM matches that passed the gathering threshold
    """

    passing_hits = []

    for hmm in current_hits.keys():
        allpass = False  # all the domains pass if the sum of bit scores meets the threshold

        if sum([bitscore for _, bitscore in current_hits[hmm].values()]) >= gathering_thresholds[hmm][0]:
            allpass = True

        for startpos in sorted(current_hits[hmm].keys()):

            matchstates, bitscore = current_hits[hmm][startpos]

            if allpass or bitscore >= gathering_thresholds[hmm][1]:
                passing_hits.append((hmm, matchstates))

    return passing_hits


####################################################################################################

def create_domsbyprot(fasta_infile,
                      path_to_pfam=REPO_DIR + 'pfam/hmms-v32/',
                      pfam_version='32',
                      results_directory=REPO_DIR + 'domains/processed-v32/',
                      concatenated_file=REPO_DIR + 'domains/allhmmresbyprot-v32.tsv.gz',
                      filtered_outfile=REPO_DIR + 'domains/domsbyprot-v32.txt.gz'):
    """
    :param fasta_infile: full path to the FASTA file that HMMER was run on (for header information)
    :param path_to_pfam: full path to a directory containing all Pfam HMMs (required for filtering results)
    :param pfam_version: version of the Pfam database we are using (default is 32)
    :param results_directory: full path to directory where tab-delimited HMMER results are stored
    :param concatenated_file: full path to a file generated by create_allhmmresbyprot()
    :param filtered_outfile: map the HMM match states to the sequence indices and residues at those indices
                    across all proteins for which we have results; note that these results are
                    guaranteed to be complete, pass default gathering thresholds, and have the
                    required amino acid at match states with high information content
    :return: None, but print success message once filtered domain results have been written to the specified
             output file
    """

    # Find all HMMs to consider:
    allhmms = find_domains_from_file(concatenated_file)

    # Determine the set of "required" states for all HMMs, where relevant
    required_states = {}  # hmm_id -> matchstate -> amino acid required
    for hmm in allhmms:
        current_req_states = get_high_information_content(path_to_pfam + hmm + '.hmm')
        if len(current_req_states.keys()) > 0:
            required_states[hmm] = current_req_states

    # Determine the lengths and instance/sequence gathering threshold cutoffs:
    hmm_lengths = {}
    gacutoff = {}
    for hmm in allhmms:
        hmm_handle = open(path_to_pfam + hmm + '.hmm')
        for hmm_line in hmm_handle:
            if hmm_line.startswith('LENG'):
                hmm_lengths[hmm] = hmm_line[:-1].split()[1]
            if hmm_line.startswith('GA'):
                gacutoff[hmm] = (float(hmm_line[:-1].split()[1].replace(';', '')),
                                 float(hmm_line[:-1].split()[2].replace(';', '')))
                break
        hmm_handle.close()

    prot_to_domlist = {}  # prot_id -> list of domains (hmm_id, set([(matchstate,index),...])) in that protein
    current_protid = ''  # current protein being processed
    current_sum = {}  # pfamID -> sum

    # Process all domains!
    res_handle = gzip.open(concatenated_file) if concatenated_file.endswith('.gz') else open(concatenated_file)
    for res_line in res_handle:
        if res_line.startswith('#') or len(res_line.strip().split('\t')) < 10:
            continue

        (prot_id, hmm_id, evalue, bit_score, targ_start, targ_end, hmm_seq,
         targ_seq, hmm_pos, desc) = res_line[:-1].split('\t')[:10]

        # Save all relevant domain hits for the previous protein:
        if prot_id != current_protid:
            passing_hits = return_passing_hits(current_sum, gacutoff)

            if len(passing_hits) > 0:
                if current_protid not in prot_to_domlist:
                    prot_to_domlist[current_protid] = []
                for hit in passing_hits:
                    prot_to_domlist[current_protid].append(hit)

            current_protid = prot_id
            current_sum = {}

        # (1) Include only COMPLETE domains:
        if not hmm_pos.startswith('1,') or not hmm_pos.endswith(',' + hmm_lengths[hmm_id]):
            continue

        # (2) Make sure that the first and last positions are ungapped:
        mstate_to_seq = zip(hmm_pos.split(','), list(targ_seq))
        if mstate_to_seq[-1][1] == '-' or mstate_to_seq[0][1] == '-':
            continue

        # (3) High information-content sites must have the appropriate residue assignment:
        if hmm_id in required_states:
            bad_match = False
            for hmm_state, seq_aa in mstate_to_seq:
                if hmm_state in required_states[hmm_id] and seq_aa != required_states[hmm_id][hmm_state]:
                    bad_match = True
                    break
            if bad_match:
                continue

        # Map HMM match state -> the HMM sequence expected -> the actual sequence hit
        seq = range(int(targ_start) - 1, int(targ_end))  # HMM sequence
        i = 0  # keep track of ungapped matches
        curr_hmm_states = []
        for hmm_state, seq_aa in mstate_to_seq:
            if seq_aa != '-':
                curr_hmm_states.append((hmm_state, seq[i], seq_aa.upper()))  # THIS HAS BEEN CHECKED!! REAL (0) INDICES!
                i += 1

        # Store this domain hit (along with match state mapping) to potentially include later on.
        if hmm_id not in current_sum:
            current_sum[hmm_id] = {}
        start_pos = min([seq_index for _, seq_index, _ in curr_hmm_states])
        if start_pos not in current_sum[hmm_id] or float(bit_score) > current_sum[hmm_id][start_pos][1]:
            current_sum[hmm_id][start_pos] = (tuple(curr_hmm_states), float(bit_score))

    res_handle.close()

    # Store information for the very last protein:
    if len(current_sum.keys()) > 0:
        passing_hits = return_passing_hits(current_sum, gacutoff)

        if len(passing_hits) > 0:
            if current_protid not in prot_to_domlist:
                prot_to_domlist[current_protid] = []
            for hit in passing_hits:
                prot_to_domlist[current_protid].append(hit)

    # Write out all sorted results:
    output_handle = gzip.open(filtered_outfile, 'w') if filtered_outfile.endswith('gz') else open(filtered_outfile, 'w')
    output_handle.write('# All COMPLETE Pfam domains (version ' + pfam_version + '.0, ' +
                        HMMINFO[pfam_version]['date'] + ', ' + HMMINFO[pfam_version]['entries'] + ' entries) that ' +
                        'passed the gathering threshold found in all amino acid sequences from\n')
    output_handle.write('# ' + fasta_infile + '\n')
    output_handle.write('# Original, unfiltered by-domain HMMER results found in ' + results_directory + '\n')
    output_handle.write('\t'.join(['#Protein_Sequence_ID', 'Pfam_HMM_ID', 'matchstate:AA-0-index:AA-value']) + '\n')

    final_domain_count = 0
    for prot_id in sorted(prot_to_domlist.keys()):
        for m in sorted(prot_to_domlist[prot_id]):
            output_handle.write(prot_id + '\t' + m[0] + '\t' +
                                ','.join([str(a[0]) + ':' + str(a[1]) + ':' + str(a[2]) for a in m[1]]) + '\n')
            final_domain_count += 1

    output_handle.close()
    sys.stderr.write('Condensed ' + concatenated_file + ' into ' + filtered_outfile + '\n')
    sys.stderr.write(str(final_domain_count) + ' total (1) complete, (2) non-deprecated domains that ' +
                     '(3) passed the gathering threshold!\n')


####################################################################################################

if __name__ == "__main__":

    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Concatenate and filter all processed HMMER domain hit results.')

    parser.add_argument('--pfam_path', type=str, default=REPO_DIR + 'pfam/',
                        help='Full path to a directory where Pfam HMMs should be stored')
    parser.add_argument('--pfam_version', type=str, default='32', choices=[str(n) for n in range(28, 33)],
                        help='Pfam version we are running on.')

    parser.add_argument('--fasta_infile', type=str, help='Full path to fasta-formatted sequence file to run HMMER on.',
                        default=REPO_DIR + 'human_test_sequences.fa')
    parser.add_argument('--results_path', type=str, default=REPO_DIR + 'domains/',
                        help='Full path to a directory where domain search results will be stored')
    parser.add_argument('--hmmer_results', type=str, default=REPO_DIR + 'domains/all-hmmer-results-by-prot-v32.txt.gz')
    parser.add_argument('--processed_results', type=str, default=REPO_DIR + 'domains/all-domains-by-prot-v32.txt.gz')

    parser.add_argument('--concatenate_hmmer_results', dest='concatenate_hmmer_results', action='store_true',
                        default=False,
                        help='Concatenate individual HMMER results file into one, nonredundant results file.')
    parser.add_argument('--filter_domains', dest='filter_domains', action='store_true', default=False,
                        help='Filter domains from concatenated results file by domain length and quality.')

    args = parser.parse_args()

    # edit path names if need be
    if not args.pfam_path.endswith('/'):
        args.pfam_path = args.pfam_path + '/'
    if not args.results_path.endswith('/'):
        args.results_path = args.results_path + '/'

    # automatically set output file names based on input arguments:
    concatenated_results_file = args.hmmer_results
    filtered_results_file = args.processed_results

    if not args.concatenate_hmmer_results and not args.filter_domains:
        sys.stderr.write(
            'No function specified, please run:\n' +
            ' '.join(['python', 'create_domain_output.py',
                      '--concatenate_hmmer_results', '--filter_domains',
                      '--pfam_path', args.pfam_path,
                      '--pfam_version', args.pfam_version,
                      '--fasta_infile', args.fasta_infile,
                      '--results_path', args.results_path,
                      '--hmmer_results', concatenated_results_file,
                      '--processed_results', filtered_results_file]) + '\n'
        )
        sys.exit(1)

    if args.concatenate_hmmer_results:
        # Concatenate all results from multiple files, removing duplicates as well as we can
        create_allhmmresbyprot(args.fasta_infile,
                               args.pfam_version,
                               args.results_path,
                               concatenated_results_file)

        sys.stderr.write(
            'For final, filtered output, remember to run:\n' +
            ' '.join(['python', 'create_domain_output.py',
                      '--filter_domains',
                      '--pfam_path', args.pfam_path,
                      '--pfam_version', args.pfam_version,
                      '--fasta_infile', args.fasta_infile,
                      '--results_path', args.results_path,
                      '--hmmer_results', concatenated_results_file,
                      '--processed_results', filtered_results_file]) + '\n'
        )

    if args.filter_domains:
        # Restrict to domains that:
        # (1) are complete (i.e., matched from the very start to the very end of the HMM)
        # (2) passed the gathering threshold (taking into account both domain- and sequence-based cutoffs)
        # (3) have the appropriate residue at high information content positions (to remove "deprecated" domains)
        if not os.path.isfile(concatenated_results_file):
            sys.stderr.write(
                'No such file: ' + concatenated_results_file + '! Please run:\n' +
                ' '.join(['python', 'create_domain_output.py',
                          '--concatenate_hmmer_results',
                          '--pfam_path', args.pfam_path,
                          '--pfam_version', args.pfam_version,
                          '--fasta_infile', args.fasta_infile,
                          '--results_path', args.results_path,
                          '--hmmer_results', concatenated_results_file]) + '\n'
            )

        create_domsbyprot(args.fasta_infile,
                          args.pfam_path + 'hmms-v' + args.pfam_version + '/',
                          args.pfam_version,
                          args.results_path,
                          concatenated_results_file,
                          filtered_results_file)
