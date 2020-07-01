#!/usr/bin/python

"""
Download latest version of Pfam and set up standard directories to run processhmmer.py
Contact snadimpa@princeton.edu with questions.
"""

import os
import gzip
import sys
import argparse
from subprocess import call


####################################################################################################

def pfam_release_info(current_release_loc='ftp://ftp.ebi.ac.uk/pub/databases/Pfam/' +
                                          'current_release/relnotes.txt'):
    """
    :param current_release_loc: complete URL to Pfam's current release notes (.txt format)
    :return: current Pfam release (e.g., 31.0), month and year of release (e.g., March 2017), and the total
             number of entries (e.g., 16,712)
    """

    # Download the current release notes and save to a temporary directory
    call(['wget', current_release_loc, '-O', 'pfam-current-release.txt'])

    rel_date, num_entries = '', ''
    with open('pfam-current-release.txt') as curr_rel_handle:
        curr_rel_handle.next()
        release = curr_rel_handle.next().strip().split()[1]  # RELEASE information always on the second line

        for rel_line in curr_rel_handle:
            if rel_line.strip().startswith(release):
                rel_date = rel_line.strip().split()[1]
                num_entries = rel_line.strip().split()[2]
                break

    # Remove the temporary file:
    call(['rm', 'pfam-current-release.txt'])

    month_abbreviations = {'01': 'January', '02': 'February', '03': 'March', '04': 'April', '05': 'May',
                           '06': 'June', '07': 'July', '08': 'August', '09': 'September', '10': 'October',
                           '11': 'November', '12': 'December'}

    month = month_abbreviations[rel_date.split('/')[0]]
    year = ('20' if int(rel_date.split('/')[-1]) < 40 else '19') + rel_date.split('/')[-1]

    return release, month + ' ' + year, "{:,}".format(int(num_entries))


####################################################################################################

def pfam_hmm_download(output_directory='hmms-v32/', pfam_version='32',
                      hmms_link='ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz'):
    """
    :param output_directory: full path to a directory where all separate HMM files will be written to
    :param pfam_version: version of Pfam that we are using (default is 32)
    :param hmms_link: complete URL to Pfam's current set of "A" HMMs
    :return: None. Download the Pfam-A HMM file, parse each individual HMM into a separate file.
    """

    # Create output directory if need be:
    if not output_directory.endswith('/'):
        output_directory = output_directory + '/'
    if not os.path.isdir(output_directory):
        call(['mkdir', output_directory])

    # Download the file containing all Pfam-A HMMs of interest:
    temporary_hmm_file = 'Pfam-A.current-v' + pfam_version + '-release.hmm' + (
        '.gz' if hmms_link.endswith('gz') else '.txt')
    call(['wget', hmms_link, '-O', temporary_hmm_file])

    # Keep track of current HMM accession and name (to print out at the appropriate time)
    current_name, current_accession = '', ''

    # Keep track of all lines belonging to this HMM:
    current_lines = []

    hmm_infile = gzip.open(temporary_hmm_file) if temporary_hmm_file.endswith('gz') else open(temporary_hmm_file)

    for l in hmm_infile:
        current_lines.append(l)

        if l.strip() == '//':
            if len(current_lines) > 0:
                hmm_outfile = open(output_directory + current_accession.split('.')[0] + '_' + current_name + '.hmm',
                                   'w')
                map(hmm_outfile.write, current_lines)
                hmm_outfile.close()
            current_name, current_accession = '', ''
            current_lines = []

        else:
            if l.startswith('NAME '):
                current_name = l.strip().split()[-1]
            if l.startswith('ACC '):
                current_accession = l.strip().split()[-1]
    hmm_infile.close()

    call(['rm', temporary_hmm_file])

    sys.stderr.write('Finished processing all HMMs from ' + hmms_link + ' into ' + output_directory + '\n')


####################################################################################################

if __name__ == "__main__":

    # parse the command-line arguments
    script_path = os.path.dirname(os.path.abspath(__file__)) + '/'

    parser = argparse.ArgumentParser(description='Download all HMMs from the most recent version of Pfam.')
    parser.add_argument('--pfam_path', type=str, help='Full path to a directory where Pfam HMMs should be stored.',
                        default=script_path + 'pfam/')
    args = parser.parse_args()

    if not args.pfam_path.endswith('/'):
        args.pfam_path = args.pfam_path + '/'

    # get information about the most recent version of Pfam:
    version, date, entries = pfam_release_info()
    current_pfam_version = version[:version.rfind('.')]
    sys.stderr.write('Pfam version ' + version + ', released ' + date + ', contains ' + entries + ' total HMMs.\n')

    for directory in [args.pfam_path,
                      args.pfam_path + 'hmms-v' + current_pfam_version]:
        if not os.path.isdir(directory):
            call(['mkdir', directory])

    # Download the most recent version of Pfam, storing HMMs in the proper directory:
    pfam_hmm_download(args.pfam_path + 'hmms-v' + current_pfam_version + '/', current_pfam_version)
