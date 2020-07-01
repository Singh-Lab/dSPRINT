#!/usr/bin/python

"""
Call relevant functions from the hmmer.py script in order to find matches to 
particular Pfam domains (and process the output from HMMER).

Contact snadimpa@princeton.edu with questions.
"""

import os
import sys
import gzip
import argparse
import time
from subprocess import call
from string import ascii_letters, digits
from random import choice
from Bio import pairwise2
from hmmer import find_domains


####################################################################################################

def idgen(size=20, chars=ascii_letters + digits):
    """
    :param size: length of string to generate
    :param chars: characters to randomly draw from in resulting string
    :return: string of length size containing random characters drawn from chars with replacement
             (intended use as a temporary file name); lowercase & uppercase letters and digits 0-9
    """

    return ''.join(choice(chars) for _ in xrange(size))


####################################################################################################

def process_hmmer_output(hmmresfile, outputfile, error_logfile, pfam):
    """
    :param hmmresfile: full path to the results file from running finddom from hmmer.py
    :param outputfile: full path to the output file containing a list of complete matches
    :param error_logfile: full path to the output file to log runtime error information
    :param pfam: dictionary of Pfam HMM ID -> the representative sequence for that HMM (obtained from the
                 raw .hmm files)
    :return: None. Output formatted results to the specified outfile. This entire process is necessary
             to determine WHICH match states aligned to which sequence positions (older version of HMMER
             makes this non-trivial)
    """

    output_handle = gzip.open(outputfile, 'w') if outputfile.endswith('gz') else open(outputfile, 'w')

    input_handle = gzip.open(hmmresfile) if hmmresfile.endswith('gz') else open(hmmresfile)

    # Skip the infile header and write out a new header:
    input_handle.next()
    output_handle.write('\t'.join(['#Target', 'HMM', 'E-value', 'BitScore', 'TargetStart', 'TargetEnd', 'HMM_Seq',
                                   'Ali_Seq', 'Pos_Seq', 'Description']) + '\n')

    for current_line in input_handle:
        # Tab-delineated file generated from the "finddom" function in hmmer.py:
        prot_id, hmm_id, evalue, bitscore, tstart, tend, hmmmatch, alimatch, desc = current_line[:-1].split('\t')[:9]

        if hmm_id not in pfam:
            sys.stderr.write('Error: Could not find ' + hmm_id + ' in pfam dictionary!\n')
            continue

        # Reformat the hmmmatch sequence off the bat; usually the "+" means the match was positive at
        #  that position in the HMM alignment (says nothing about the complete match). Full description
        #  of HMMER output can be found in the User's Guide:
        # http://eddylab.org/software/hmmer3/3.1b2/Userguide.pdf
        hmmmatch = hmmmatch.replace('+', '')

        # Sometimes only part of the domain matched. Figure out the starting index of the match:
        startingindex = pfam[hmm_id][0].find(hmmmatch.replace('.', '').replace('+', ''))

        # If the resulting HMM match doesn't match what we expect for the overall Pfam sequence match, then
        #  we will not be able to find the starting index. In this case, we need to "translate" our HMM match
        #  into what we expected to see via a quick alignment:
        if startingindex < 0 or len(hmmmatch) != len(alimatch):

            # Let's align our hmmmatch to what we SHOULD have, the Pfam sequence.
            # A mismatch is only -0.5 whereas gaps are -2, so mismatch is preferred.
            newh, alim = pairwise2.align.globalms(hmmmatch.replace('.', ''), pfam[hmm_id][0], 2, -0.5, -1.5, -1)[0][:2]

            newhmmmatch = []
            i, j, k = 0, 0, 0  # indices into original HMM match, aligned HMM match, and aligned expected Pfam sequence

            # ....Simple quick alignment adjustment:
            while i < len(hmmmatch) or j < len(newh) or k < len(alim):  # alim must ALWAYS be >= in length.
                if newh[j] == '-' and alim[k] != '-' and i >= len(hmmmatch):  # gap at the END, add nothing & move on
                    j += 1
                    k += 1
                elif hmmmatch[i] == '.':
                    newhmmmatch.append('.')
                    i += 1
                elif newh[j] == '-' and alim[k] != '-' and i == 0:  # gap at the BEGINNING, add nothing & move on
                    j += 1
                    k += 1
                elif newh[j] == '-' and alim[k] != '-' and i < len(hmmmatch):  # ADD THIS CHARACTER
                    newhmmmatch.append(alim[k])
                    j += 1
                    k += 1
                elif newh[j] != '-' and newh[j] == hmmmatch[i] and alim[k] == '-':  # DELETE THIS CHARACTER
                    i += 1
                    j += 1
                    k += 1
                elif newh[j] != '-' and newh[j] == hmmmatch[i] and alim[k] != '-':  # REPLACE this character
                    newhmmmatch.append(alim[k])
                    i += 1
                    j += 1
                    k += 1
                elif newh[j] == alim[k] == hmmmatch[i]:
                    newhmmmatch.append(hmmmatch[i])
                    i += 1
                    j += 1
                    k += 1
                else:  # No idea what this case is! Let's check it out...
                    sys.stderr.write(str(i) + ':' + hmmmatch[i] + '\t' + hmmmatch + '\n')  # Original HMM match
                    sys.stderr.write(str(j) + ':' + newh[j] + '\t' + newh + '\n')  # Aligned HMM match
                    sys.stderr.write(str(k) + ':' + alim[k] + '\t' + alim + '\n')  # Aligned Pfam sequence match

            # Did this fix the problem!?
            newhmmmatch = ''.join(newhmmmatch)
            startingindex = pfam[hmm_id][0].find(newhmmmatch.replace('.', '').replace('+', ''))

            if startingindex < 0 or len(newhmmmatch) != len(alimatch):
                logfile_handle = open(error_logfile, 'a')
                logfile_handle.write(hmm_id + '\n')
                logfile_handle.write('Orig HMM Match:\t' + hmmmatch + '\n')
                logfile_handle.write('Align HMM Matc:\t' + newh + '\n')
                logfile_handle.write('Align Pfam Mat:\t' + alim + '\n')
                logfile_handle.write('New HMM Match: \t' + newhmmmatch + '\n')
                logfile_handle.write('Orig Ali Match:\t' + alimatch + '\n\n')
                logfile_handle.close()
                sys.stderr.write('Problem in ' + hmm_id + '\n')
                continue

            else:
                hmmmatch = newhmmmatch

        # Actual mapping of Pfam indices:
        mapindices = {i + 1: int(k) for i, k in enumerate(pfam[hmm_id][1].split(','))}

        startingindex += 1  # Results should be 1-indexed, not 0-indexed, for consistency with HMMER.

        findex = []
        currgap = True
        currgapi = 0

        # Determine the actual match state -> sequence mapping. Insertion states are prefixed with an "a"
        for i in xrange(len(hmmmatch)):
            if hmmmatch[i] != '.':
                currgap = False
                findex.append(str(mapindices[startingindex]))
                startingindex += 1
            elif hmmmatch[i] == '.':
                if currgap:
                    currgapi += 1
                else:
                    currgap = True
                    currgapi = 0
                findex.append('a' + str(mapindices[startingindex] - 1) + '-' + str(currgapi))

        output_handle.write('\t'.join([prot_id, hmm_id, evalue, bitscore, tstart, tend,
                                       hmmmatch, alimatch, ','.join(findex), desc]) + '\n')

    input_handle.close()
    output_handle.close()


####################################################################################################

def get_pfamseq(hmm_file):
    """
    :param hmm_file: full path to a Pfam HMM file
    :return: consensus sequence representing the HMM, needed to parse hmmsearch results (v2.3.2) to
             get the match state -> sequence index -> amino acid mapping
    """

    pfamseq = []
    for current_line in open(hmm_file):
        i = current_line.strip().split()
        if len(i) == 26 and i[21].isdigit():
            pfamseq.append(i[22])
    pfamseq = ''.join(pfamseq)

    return pfamseq


####################################################################################################

def find_domain_matches(hmm_file, output_file, infile, error_logfile, ids=()):
    """
    :param hmm_file: full path to an HMM. Remember our naming convention MUST be PATH + PfamID_PfamName.hmm
    :param output_file: full path to an output file to store the output
    :param infile: full path to a FASTA file with protein sequences that we want to find domains in
    :param error_logfile: full path to the output file to log runtime error information
    :param ids: a subset of sequence identifiers that we care about from infile
    :return: find all sequence matches to the HMM profile in hmmfile for the specified sequence IDs found in
             the input fasta file. Write the *complete* (not necessarily all confident), easily-parsed results
             to the output file.
    """

    # Pfam "consensus" sequence is needed to parse hmmsearch results; obtained from the raw .hmm file:
    pfamseq = get_pfamseq(hmm_file)

    chmm = hmm_file.split('/')[-1].replace('.hmm', '')  # PfamID_PfamName

    output_handle = gzip.open(output_file, 'w') if output_file.endswith('gz') else open(output_file, 'w')
    output_handle.write('# All matching hits from the HMM found in ' + hmm_file + '\n')
    output_handle.write('# on the protein sequences found in ' + infile + '\n')
    output_handle.write(
        '\t'.join(['#TargetID', 'HMM_Name', 'E-value', 'BitScore', 'TargetStart', 'TargetEnd', 'HMM_Seq',
                   'Target_Seq', 'HMM_Pos', 'Description']) + '\n')

    # NOW, let's search, seq by seq...
    input_handle = gzip.open(infile) if infile.endswith('gz') else open(infile)
    while True:
        current_line = input_handle.readline()
        if not current_line:
            break

        # Check if this is a sequence we actually want to run on:
        if current_line.startswith('>'):
            sequence_id = current_line[1:-1].split()[0]
            if len(ids) > 0 and sequence_id not in ids:
                continue

            tmpfiles = ['/tmp/' + idgen() + '.seq',  # sequence input file (single sequence)
                        '/tmp/' + idgen() + '.hmmres',  # original unparsed output
                        '/tmp/' + idgen() + '.hmmres1']  # complete parsed output

            # Write out current header and following sequence as input for hmmsearch
            current_line_position = input_handle.tell()  # keep track of our current file location
            next_line = input_handle.readline()  # and read in the following line
            if not next_line:
                input_handle.seek(current_line_position)  # go back one, and reset.
                continue
            sequence_handle = open(tmpfiles[0], 'w')
            sequence_handle.write(current_line)

            while not next_line.startswith('>'):
                sequence_handle.write(next_line.strip())
                current_line_position = input_handle.tell()
                next_line = input_handle.readline()
                if not next_line:
                    break
            input_handle.seek(current_line_position)  # reset to the last line for beginning of loop.
            sequence_handle.write('\n')
            sequence_handle.close()

            # Now, RUN hmmsearch for EACH hmmID in our list.
            find_domains([hmm_file], tmpfiles[0], tmpfiles[1])

            os.system('rm ' + tmpfiles[0])

            # CHECK IF THIS FAILED!!!
            if sum(1 for _ in open(tmpfiles[1])) < 2:
                sys.stderr.write('Failed to find any HMM matches for ' + hmm_file + ' in ' + sequence_id + '.\n')
                os.system('rm ' + tmpfiles[1])
                continue

            # we need the "pfam sequence" for our particular HMM so that we can match positions to match states.
            process_hmmer_output(tmpfiles[1], tmpfiles[2], error_logfile,
                                 {chmm[chmm.find('_') + 1:]: (pfamseq, ','.join(map(str, range(1, len(pfamseq) + 1))))})

            os.system('rm ' + tmpfiles[1])

            with open(tmpfiles[2]) as y:
                y.next()
                for l2 in y:
                    # Write out the exact same line, but update the HMM name to include Pfam ID (e.g., PF00096) also
                    output_handle.write(
                        l2[:-1].split('\t')[0] + '\t' + chmm + '\t' + '\t'.join(l2[:-1].split('\t')[2:]) + '\n')
            os.system('rm ' + tmpfiles[2])

    input_handle.close()
    output_handle.close()


####################################################################################################

if __name__ == "__main__":

    # parse command-line arguments:
    parser = argparse.ArgumentParser(description='Run, parse, and process results from HMMER 3.0 and HMMER 2.3.2.')

    # to parallelize, we can specify a *subset* of these domains to run HMMER on.
    script_path = os.path.dirname(os.path.abspath(__file__)) + '/'

    parser.add_argument('--start', type=int, default=0,
                        help='Starting 0-index of subset of domains to run on.')
    parser.add_argument('--end', type=int,
                        help='Ending 0-index of subset of domains to run on.')

    parser.add_argument('--pfam_path', type=str, default=script_path + 'pfam/',
                        help='Full path to a directory where Pfam HMMs should be stored')
    parser.add_argument('--pfam_version', type=str, default='32', choices=[str(n) for n in range(28, 33)],
                        help='Pfam version we are running on.')
    parser.add_argument('--fasta_infile', type=str, help='Full path to fasta-formatted sequence file to run HMMER on.')
    parser.add_argument('--results_path', type=str, default=script_path + 'domains/',
                        help='Full path to a directory where domain search results will be stored')

    args = parser.parse_args()

    # edit path names
    if not args.pfam_path.endswith('/'):
        args.pfam_path = args.pfam_path + '/'
    if not args.results_path.endswith('/'):
        args.results_path = args.results_path + '/'

    # (1) create output directories:
    for directory in [args.results_path,
                      # unprocessed HMMER results will go here (keep in case of a crash / for debugging):
                      args.results_path + 'hmmres-v' + args.pfam_version,
                      # final, processed output will go here, for each domain (e.g., PF00096_zf-C2H2-v32.hmmres.gz):
                      args.results_path + 'processed-v' + args.pfam_version]:
        if not os.path.isdir(directory):
            exit_code = call(['mkdir', directory])
            if exit_code != 0:
                sys.stderr.write('Could not create directory: ' + directory + '\n')
                sys.exit(1)

    # (2) find all possible hmms to run on:
    #     NOTE: we assume files are named as PfamID_PfamName.hmm (e.g., PF00096_zf-C2H2.hmm)
    hmms = sorted([a.replace('.hmm', '') for a in os.listdir(args.pfam_path + 'hmms-v' + args.pfam_version)
                   if a.startswith('PF') and a.endswith('.hmm')])
    total_hmms_available = len(hmms)

    # (3) customize log file name if need be and subset to appropriate range of HMMs
    #     NOTE: in case we had a match that we couldn't understand nor rescue, keep track of it in the logfile
    logfile = args.results_path + 'processed-v' + args.pfam_version + '/problems.log'
    if args.start != 0 or args.end < total_hmms_available:
        logfile = logfile.replace('.log', '-' + str(args.start) + '-' + str(args.end) + '.log')
        hmms = hmms[args.start:min(args.end, total_hmms_available)]  # Subset to a range of HMMs if it was specified
    fasta_handle = gzip.open(logfile, 'w') if logfile.endswith('gz') else open(logfile, 'w')
    fasta_handle.close()  # clear the logfile to keep track of only the most recent results:

    # (4) read in sequence identifiers from input FASTA file
    #     NOTE: hmmsearch causes a fuss if this file is zipped at all, so make sure it is not.
    if not os.path.isfile(args.fasta_infile):
        sys.stderr.write('Could not open specified infile: ' + args.fasta_infile + '\n')
        sys.exit(1)
    if args.fasta_infile.endswith('.gz'):
        call(['gzip', '-d', args.fasta_infile])

    allseqs = set()
    fasta_handle = gzip.open(args.fasta_infile) if args.fasta_infile.endswith('gz') else open(args.fasta_infile)
    for seq_line in fasta_handle:
        if seq_line.startswith('>'):
            allseqs.add(seq_line[1:-1].split()[0])
    fasta_handle.close()

    # (5) finally, for each HMM
    for hmm in hmms:
        hmmfile = args.pfam_path + 'hmms-v' + args.pfam_version + '/' + hmm + '.hmm'

        # (5a) check if HMM exists
        if not os.path.isfile(hmmfile):
            sys.stderr.write('\n'.join(['Could not find HMM file: ' + hmmfile,
                                        'Download latest version from Pfam:',
                                        'wget -O ' + hmmfile + ' http://pfam.xfam.org/family/' +
                                        hmm.split('_')[0] + '/hmm']) + '\n')
            continue

        # (5b) run hmmsearch to subset the genes we want to actually find domain hits in
        #     NOTE: this is just an efficiency step, as this process is FAST but gives us less information
        call(' '.join(['hmmsearch',
                       '-o /dev/null',
                       '--domtblout ' + args.results_path + 'hmmres-v' + args.pfam_version + '/' + hmm + '.hmmres-orig',
                       '-T 0',
                       '--domT 0',
                       '--incT 0',
                       '--incdomT 0',  # No cutoffs guarantees more thorough hits
                       hmmfile,
                       args.fasta_infile]), shell=True)

        # (5c) Set the subset of sequence IDs to look through (to speed up process of finding matches):
        whichseqs = set()
        if os.path.isfile(args.results_path + 'hmmres-v' + args.pfam_version + '/' + hmm + '.hmmres-orig'):
            for seq_line in open(args.results_path + 'hmmres-v' + args.pfam_version + '/' + hmm + '.hmmres-orig'):
                if seq_line.startswith('#'):
                    continue
                whichseqs.add(seq_line.strip().split()[0])

        # (5d) run on all sequences instead of subset if no hits were found (debugging purposes ONLY)
        if False and len(whichseqs) < 1:
            sys.stderr.write('Failed to get preliminary results for ' + hmm + '! ' +
                             'Running on all ' + str(len(allseqs)) + ' sequences...\n')
            whichseqs = allseqs

        # (5e) find the full matches (match state -> sequence index -> sequence residue information):
        if len(whichseqs) > 0:
            sys.stderr.write('Processing ' + hmm + '-v' + args.pfam_version + '...')

            # start the clock to measure performance
            start = time.time()

            outfile = args.results_path+'processed-v'+args.pfam_version+'/'+hmm+'-v'+args.pfam_version+'.hmmres.gz'
            find_domain_matches(hmmfile, outfile, args.fasta_infile, logfile, whichseqs)

            # end the clock and print total elapsed time:
            end = time.time()
            total_seconds = end - start
            m, s = divmod(total_seconds, 60)
            h, m = divmod(m, 60)
            d, h = divmod(h, 24)
            sys.stderr.write(hmm + ': processed ' + str(len(whichseqs)) + ' sequences in ' +
                             ':'.join(map(lambda x: str(int(x)).zfill(2), [d, h, m, s])) + '!\n')
