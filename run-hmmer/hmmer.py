#!/usr/bin/python

"""
Run HMMER 2.3.2 and/or HMMER 3 and parse output into tab-delineated files.

Archived versions of HMMER 2.3.2 and HMMER 3.0 are available for download from
  http://hmmer.janelia.org/software/archive

When using HMMER, you should cite
  Finn RD, Clements J, Eddy SR (2011) HMMER web server: interactive sequence
  similarity searching. Nucl Acids Res 39(suppl 2): W29-W37. 
  doi: 10.1093/nar/gkr367

The following programs should be available for command line use:*
  hmmsearch    -- HMMER 3.0 version
  hmmconvert   -- HMMER 3.0 version
  hmmsearch232 -- HMMER 2.3.2 version (and named as such)

*If these programs are not available or cannot be run as listed, instances in 
 the code where these programs are called will need to be changed.

Contact snadimpa@princeton.edu with questions.
"""

import os
import sys
import gzip
from random import choice
from string import ascii_letters, digits
from subprocess import call


####################################################################################################
# Extract information from HMM and sequence files
####################################################################################################

def hmmname_from_hmmfile(hmm_file):
    """
    :param hmm_file: full path to a Pfam HMM file
    :return: the name of the HMM parsed from the file
    """

    if not os.path.isfile(hmm_file):
        sys.stderr.write('Could not open ' + hmm_file + '\n')
        return

    hmm_handle = gzip.open(hmm_file) if hmm_file.endswith('gz') else open(hmm_file)

    hmm_name = 'N/A'
    for hmm_line in hmm_handle:
        if hmm_line.startswith('NAME'):
            hmm_name = hmm_line.strip().split()[1]
            break
    else:
        sys.stderr.write('Improperly formatted HMM file (' + hmm_file + '), could not find name.\n')

    hmm_handle.close()
    return hmm_name


####################################################################################################

def seqdesc_from_fasta(fasta_file):
    """
    :param fasta_file: full path to a fasta-formatted sequence file
    :return: a dictionary of sequence ID -> FASTA descriptions
    """

    fasta_handle = gzip.open(fasta_file) if fasta_file.endswith('gz') else open(fasta_file)
    descriptions = {fasta_line[1:-1].split()[0]: ' '.join(fasta_line.strip().split()[1:])
                    for fasta_line in fasta_handle if fasta_line.startswith('>')}
    fasta_handle.close()

    return descriptions


####################################################################################################
# Parse HMMER 2.3.2 and HMMER 3.1b2 results
####################################################################################################

def parsehmmer232(unparsed_output, hmm_name, seqid_to_desc):
    """
    :param unparsed_output: full path to unparsed output from HMMER 2.3.2
    :param hmm_name: Pfam HMM ID and name (e.g., PF00096_zf-C2H2)
    :param seqid_to_desc: dictionary of sequence ID -> sequence description from the original FASTA file
    :return: a tab-delimited table with all domain hits and additional important information
    """

    tmpfile = '/tmp/' + ''.join(choice(ascii_letters + digits) for _ in xrange(10)) + '.txt'
    output_handle = open(tmpfile, 'w')

    # Some variables to be used while looping through file
    reach = False  # Have we reached a match yet?
    multihmm = False  # Does the HMM match span multiple lines?
    cid = ''  # Keep track of the current ID once we've reached a match

    # Set blank variables of all fields that will be updated properly if files are properly formatted:
    prot_id, bitscore, evalue, alistart, aliend, alimatch = '', '', '', '', '', ''
    hmmmatch, hmmstart, hmmend = '', '', ''

    hmm_position = {}

    unparsed_handle = open(unparsed_output)
    for output_line in unparsed_handle:
        i = output_line.strip().split()
        if len(i) < 1:  # No pertinent information here
            continue

        # We really need the HMM start and end positions also, though...
        if len(i) > 2 and i[0] in seqid_to_desc and i[1].count('/') == 1 and \
           i[1].split('/')[0].isdigit() and i[1].split('/')[1].isdigit():
            a1 = max(len(i[0]), len('Sequence')) + 1
            hmm_position[i[1]] = (output_line[a1 + 23:a1 + 28].strip(), output_line[a1 + 29:a1 + 34].strip())

            if not hmm_position[i[1]][0].isdigit() or not hmm_position[i[1]][1].isdigit():
                sys.stderr.write(
                    '\n'.join([unparsed_output, output_line, hmm_position[i[1]][0], hmm_position[i[1]][1]]) + '\n')
                sys.exit(1)

        # If i[0] is an ID that we're interested in:
        elif len(i) > 2 and i[0][:-1] in seqid_to_desc and i[0][-1] == ':' and i[1] == 'domain':
            assert not reach, 'Could not correctly parse ' + unparsed_output

            reach = True  # We have reached a match

            # This bit is strange!! IF your ID is less than 10 characters, a ':' is tacked
            # on the end of it. If your ID is longer than 10 characters, it is truncated to
            # 10 characters WHEN USED LATER IN THE ALIGNMENT MATCH (here, it is always complete)

            currkey = i[2] + '/' + i[4].replace(',', '')
            hmmstart, hmmend = hmm_position[currkey] if currkey in hmm_position else ('', '')

            cid = i[0][:10]
            if cid.endswith(':') and len(cid) == len(i[0]):
                cid = cid[:-1]

            # Parse all possible info from this line:
            prot_id = i[0][:-1]  # Remove the ending ':'
            bitscore = i[10][:-1]  # Remove the ending ':'
            evalue = i[13]
            alistart = i[6]
            aliend = i[8][:-1]  # Remove the ending ':'
            alimatch, hmmmatch = '', ''

        elif reach:
            if i[0][:3] == '*->':  # Start of HMM match
                if '<' not in i[0]:  # Multi-line HMM match
                    multihmm = True
                    hmmmatch = i[0][3:]
                else:
                    hmmmatch = i[0][3:i[0].find('<')]

            elif i[0] == cid and len(i) > 3:  # Part of an alignment match
                if len(alimatch) < 1:
                    alimatch = i[2]
                else:
                    alimatch += i[2]  # (multi-line)

                if aliend == i[3]:  # Reached the end
                    output_handle.write('\t'.join([prot_id, hmm_name, evalue, bitscore, alistart, aliend, hmmmatch,
                                                   alimatch,
                                                   (seqid_to_desc[prot_id] + ' ' if prot_id in seqid_to_desc else '') +
                                                   'HMMStart=' + str(hmmstart) + '; HMMEnd=' + str(
                                                       hmmend) + ';']) + '\n')
                    reach = False
                    multihmm = False

            elif multihmm and (len(i) == 1 or '+' not in output_line):  # Last line of a multi-line match?
                if '<' in i[0]:
                    multihmm = False
                    hmmmatch += i[0][:i[0].find('<')]
                elif i[0] != 'CS':  # This stands for 'consensus structure' line -- we ignore that
                    hmmmatch += i[0]
    unparsed_handle.close()
    output_handle.close()

    # Overwrite the original HMMER output:
    final_output_handle = open(unparsed_output, 'w')
    final_output_handle.write('\t'.join(['#Target', 'HMM', 'E-value', 'BitScore', 'TargetStart', 'TargetEnd',
                                         'HMM_Seq', 'Ali_Seq', 'Description']) + '\n')
    with open(tmpfile) as tmp_handle:
        for output_line in sorted([output_line for output_line in tmp_handle]):
            final_output_handle.write(output_line)
    final_output_handle.close()
    call(['rm', tmpfile])


####################################################################################################

def parsehmmer3(unparsed_output, hmm_name, seqid_to_desc):
    """
    :param unparsed_output: full path to unparsed output from HMMER 3.1b2
    :param hmm_name: Pfam HMM ID and name (e.g., PF00096_zf-C2H2)
    :param seqid_to_desc: dictionary of sequence ID -> sequence description from the original FASTA file
    :return: a tab-delimited table with all domain hits and additional important information
    """

    tmpfile = '/tmp/' + ''.join(choice(ascii_letters + digits) for _ in xrange(10)) + '.txt'
    output_handle = open(tmpfile, 'w')

    reach = False  # Have we reached a match yet?
    prot_id, bitscore, evalue = '', '', ''
    hmmstart, hmmend, hmmmatch, alistart, aliend, alimatch = '', '', '', '', '', ''

    unparsed_handle = open(unparsed_output)
    for output_line in unparsed_handle:
        i = output_line.strip().split()
        if len(i) < 1:  # No pertinent information here
            continue

        if i[0] == '==':
            # We've reached one domain already, so print it out.
            if reach:
                output_handle.write(
                    '\t'.join([prot_id, hmm_name, evalue, bitscore, alistart, aliend, hmmmatch, alimatch,
                               (seqid_to_desc[prot_id] + ' ' if prot_id in seqid_to_desc else '') +
                               'HMMStart=' + str(hmmstart) + '; HMMEnd=' + str(hmmend) + ';']) + '\n')
                hmmstart, alistart, hmmmatch, alimatch, prot_id = '', '', '', '', ''

            reach = True
            # Example line:
            # == domain 2 score: 40.7 bits; conditional E-value: 1.3e-14
            bitscore = i[4]
            evalue = i[8]

        elif reach:  # Otherwise, parse out information to augment domain
            if i[0].startswith('>>') or (i[0].startswith('Internal') and not i[0].startswith(hmm_name)):
                output_handle.write(
                    '\t'.join([prot_id, hmm_name, evalue, bitscore, alistart, aliend, hmmmatch, alimatch,
                               (seqid_to_desc[prot_id] + ' ' if prot_id in seqid_to_desc else '') +
                               'HMMStart=' + str(hmmstart) + '; HMMEnd=' + str(hmmend) + ';']) + '\n')
                reach = False
                hmmstart, alistart, hmmmatch, alimatch, prot_id = '', '', '', '', ''

            elif i[0].startswith(hmm_name) and len(i[0]) <= len(hmm_name):
                # Example line:
                # fn3 2 saPenlsvsevtstsltlsWsppkdgggpitgYeveyqekgegeewqevtvprtttsvtltgLepgteYefrVqavngagegp 84
                if hmmstart == '':
                    hmmstart = i[1]
                    hmmmatch = i[2]
                else:
                    hmmmatch += i[2]
                hmmend = i[3]  # this will be reset with multi-line HMMs

            elif i[0] in seqid_to_desc:
                # Example line:
                # 7LESS_DROME 439 SAPVIEHLMGLDDSHLAVHWHPGRFTNGPIEGYRLRLSSSEGNA-TSEQLVPAGRGSYIFSQLQAGTNYTLALSMINKQGEG 520
                if prot_id == '':
                    prot_id = i[0]
                    alistart = i[1]
                    alimatch = i[2]
                else:
                    alimatch += i[2]
                aliend = i[3]  # this will ALSO be reset with multi-line alignment matches
    unparsed_handle.close()
    output_handle.close()

    # Overwrite the original HMMER output:
    final_output_handle = open(unparsed_output, 'w')
    final_output_handle.write('\t'.join(['#Target', 'HMM', 'E-value', 'BitScore', 'TargetStart', 'TargetEnd',
                                         'HMM_Seq', 'Ali_Seq', 'Description']) + '\n')
    with open(tmpfile) as tmp_handle:
        for output_line in sorted([output_line for output_line in tmp_handle]):
            final_output_handle.write(output_line)
    final_output_handle.close()
    call(['rm', tmpfile])


####################################################################################################
# Run HMMER 2.3.2 and HMMER 3.1b2
####################################################################################################

def runhmmer232(hmm_file, input_prot_file, results_out, hmm_name='', seqid_to_desc=None):
    """
    :param hmm_file: full path to a Pfam-formatted HMM
    :param input_prot_file: full path to an *unzipped* FASTA-formatted sequence file
    :param results_out: full path to where parsed, formatted output should be stored
    :param hmm_name: Pfam HMM ID and name (e.g., PF00096_zf-C2H2)
    :param seqid_to_desc: dictionary of sequence ID -> sequence description from the original FASTA file
    :return: none, but run HMMER 2.3.2, parse the results, and store those results in results_out
    """

    # Check if all files exist:
    for infile in [hmm_file, input_prot_file]:
        if not os.path.isfile(infile):
            sys.stderr.write('Could not open file ' + infile + '.\n')
            return
    if not os.path.isdir('/'.join(results_out.split('/')[:-1])):
        sys.stderr.write('Could not create outputfile in ' +
                         '/'.join(results_out.split('/')[:-1]) + ', no such directory.\n')
        return

    # Extract hmmname and descriptions if not provided:
    if hmm_name == '':
        hmm_name = hmmname_from_hmmfile(hmm_file)
    if not seqid_to_desc:
        seqid_to_desc = seqdesc_from_fasta(input_prot_file)

    with open(hmm_file) as x:
        firstline = x.next()

    if not firstline.startswith('HMMER2.0'):
        tmpfile = '/tmp/' + ''.join(choice(ascii_letters + digits) for _ in xrange(10)) + '.txt'
        call('hmmconvert -2 ' + hmm_file + ' > ' + tmpfile, shell=True)
        call('chmod 755 ' + tmpfile, shell=True)
        hmm_file = tmpfile

    call("hmmsearch232 -Z 10 --domT 0 " + hmm_file + " " + input_prot_file + " > " + results_out, shell=True)

    if not firstline.startswith('HMMER2.0'):
        call(['rm', hmm_file])

    parsehmmer232(results_out, hmm_name, seqid_to_desc)


####################################################################################################

def runhmmer3(hmm_file, input_prot_file, results_out, hmm_name='', seqid_to_desc=None):
    """
    :param hmm_file: full path to a Pfam-formatted HMM
    :param input_prot_file: full path to an *unzipped* FASTA-formatted sequence file
    :param results_out: full path to where parsed, formatted output should be stored
    :param hmm_name: Pfam HMM ID and name (e.g., PF00096_zf-C2H2)
    :param seqid_to_desc: dictionary of sequence ID -> sequence description from the original FASTA file
    :return: none, but run HMMER 3.1b2, parse the results, and store those results in results_out
    """

    # Check if all files exist:
    for infile in [hmm_file, input_prot_file]:
        if not os.path.isfile(infile):
            sys.stderr.write('Could not open file ' + infile + '.\n')
            return
    if not os.path.isdir('/'.join(results_out.split('/')[:-1])):
        sys.stderr.write('Could not create output file in ' +
                         '/'.join(results_out.split('/')[:-1]) + ', no such directory.\n')
        return

    # Extract hmmname and descriptions if not provided:
    if hmm_name == '':
        hmm_name = hmmname_from_hmmfile(hmm_file)
    if not seqid_to_desc:
        seqid_to_desc = seqdesc_from_fasta(input_prot_file)

    call("hmmsearch -Z 10 --notextw --domT 0 " + hmm_file + " " + input_prot_file + " > " + results_out, shell=True)
    parsehmmer3(results_out, hmm_name, seqid_to_desc)


####################################################################################################
# Primary parser functionality
####################################################################################################

def find_domains(hmms_to_run, input_prot_file, results_out):
    """
    :param hmms_to_run: set of HMMs to run HMMER 2.3.2 and 3.1b2 on
    :param input_prot_file: full path to an *unzipped* FASTA-formatted sequence file (to run HMMER on)
    :param results_out: full path to an output file to store all unique output from HMMER
    :return: None, but write to results_out
    """
    """Runs HMMER 2.3.2 and HMMER 3 using the HMMs listed in hmmlist on the
    filenames found in the inputdir, and writes output files to the outputdir
    with the name origname + suffix + '.hmmres'."""

    # Check the type of the hmmlist!
    if type(hmms_to_run) not in [list, set, tuple]:
        hmms_to_run = [hmms_to_run]

    # Check to make sure files exist:
    if not os.path.isfile(input_prot_file):
        sys.stderr.write('In find_domains: Could not find input file ' + input_prot_file + '\n')
        return
    if not os.path.isdir(results_out[:results_out.rfind('/')]):
        sys.stderr.write(
            'In find_domains: Could not find output directory ' + results_out[:results_out.rfind('/')] + '\n')
        return
    for current_hmm in hmms_to_run:
        if not os.path.isfile(current_hmm):
            sys.stderr.write('In find_domains: Could not find HMM file ' + current_hmm + '\n')
            return

    # Get protein descriptions
    descs = seqdesc_from_fasta(input_prot_file)

    # Run both HMMER 2.3.2 and HMMER 3, parse the output, and store it all
    tmpfiles = []
    for hmmfile in hmms_to_run:
        hmmname = hmmname_from_hmmfile(hmmfile)

        for hmmfunc in [runhmmer232, runhmmer3]:
            tmpfiles.append('/tmp/' + ''.join(choice(ascii_letters + digits) for _ in xrange(10)) + '.tmp')
            hmmfunc(hmmfile, input_prot_file, tmpfiles[-1], hmmname, descs)

    # Keep ALL UNIQUE output, sort it, and print it into outfile
    lines = set()
    for f in tmpfiles:
        for l in open(f):
            if len(l.strip().split('\t')) < 9:
                continue
            if l.strip() != '':
                lines.add(l)
        call(['rm', f])

    output_handle = open(results_out, 'w')
    map(output_handle.write, sorted(list(lines)))
    output_handle.close()


####################################################################################################
# MAIN
# Note that this exists purely for *testing and debugging* purposes.
# Functions should be imported from this file to use separately.
# I offer no guarantees on this testing code!
####################################################################################################

if __name__ == "__main__":

    sys.stderr.write('Testing all functions in ' + sys.argv[0] + '...\n')

    data_path = '/home/snadimpa/datadb/'

    # To run the following set of HMMs:
    hmms = ['PF00096_zf-C2H2', 'PF02892_zf-BED', 'PF06220_zf-U1', 'PF09237_GAGA',
            'PF11931_DUF3449', 'PF12171_zf-C2H2_jaz', 'PF12756_zf-C2H2_2',
            'PF12874_zf-met', 'PF12907_zf-met2', 'PF13909_zf-H2C2_5', 'PF13912_zf-C2H2_6',
            'PF13913_zf-C2HC_2']

    hmms = [data_path + 'pfam/' + hmm + '.hmm' for hmm in hmms]

    # On the following PROTEIN input file:
    seqfile = data_path + 'flybase/pep/dmel-all-translation-r6.04.fasta.gz'

    # And store the HMMER results in the following output file:
    outfile = data_path + 'hmmer/dmel-r6.04_c2h2zf_hmmerresults.tsv'

    sys.stderr.write('Looking for matches using the following ' + str(len(hmms)) + ' HMM files:\n')
    for hmm in hmms:
        sys.stderr.write(hmm + '\n')
    sys.stderr.write('Searching sequences in ' + seqfile + ' for matches, writing output to ' + outfile + '\n')

    # We simply run:
    find_domains(hmms, seqfile, outfile)

    # If we wanted to choose a particular HMM and run only HMMER 2.3.2 or HMMER 3.0:

    sys.stderr.write('Running HMMER 2.3.2 using hmm ' + hmms[0] + '\n')
    runhmmer232(hmms[0], seqfile, outfile)
    sys.stderr.write('Success!\n')

    sys.stderr.write('Running HMMER 3.0 using hmm ' + hmms[0] + '\n')
    runhmmer3(hmms[0], seqfile, outfile)
    sys.stderr.write('Success!\n')
