#!/usr/bin/python

"""
Create necessary track files reflecting Pfam domain locations in human protein sequences

Contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import gzip
import argparse
from config import data_path, GENOME_BUILD
from subprocess import call

####################################################################################################
# CONSTANTS
####################################################################################################

PFAM_VERSION = '31'


####################################################################################################

def range_from_indices(index_list):
    """
    :param index_list: list or set of indices
    :return: string corresponding to ordered intervals of values
    """

    indices = sorted(list(set(map(int, index_list))))

    ranges = []
    current_range = []

    for i in indices:
        if len(current_range) < 1:
            current_range = [i]
        elif len(current_range) < 2:
            if i == current_range[0]+1:
                current_range.append(i)
            else:
                ranges.append(str(current_range[0]))
                current_range = [i]
        else:
            if i == current_range[1]+1:
                current_range[1] += 1
            else:
                ranges.append(str(current_range[0])+'-'+str(current_range[1]))
                current_range = [i]

    # add the final range
    if len(current_range) == 1:
        ranges.append(str(current_range[0]))
    elif len(current_range) == 2:
        ranges.append(str(current_range[0])+'-'+str(current_range[1]))

    return ranges


####################################################################################################

def create_domain_tracks(seq_file, domain_locs, out_trackfile):
    """
    :param seq_file: full path to a fasta-formatted file containing all human protein sequences
                     (to get lengths of protein isoforms required for tracks)
    :param domain_locs: full path to a tab-delimited file containing domain locations in human proteins
    :param out_trackfile: full path to a file to write out track information
    :return: None, but print success message upon write of two output files
    """

    # we want to create a new tab-delimited file with the following columns:
    # [0] Ensembl protein ID, [2] PfamID_PfamName, [3] protein length,
    # [4] comma-delimited list of range(s) that the domain spans

    print 'made it to function'

    # (1) process the domain locations per Ensembl ID:
    prot_to_domains = {}
    dom_handle = gzip.open(domain_locs, 'rt') if domain_locs.endswith('gz') else open(domain_locs)
    for dom_line in dom_handle:
        if dom_line.startswith('#'):
            continue
        protein_id, domain_name, match_states = dom_line[:-1].split('\t')[:3]
        if protein_id not in prot_to_domains:
            prot_to_domains[protein_id] = set()

        # 0-index positions that this domain "covers"
        domain_range = set([int(ms.split(':')[1]) for ms in match_states.split(',')])
        prot_to_domains[protein_id].add((domain_name, str(min(domain_range))+'-'+str(max(domain_range))))
    dom_handle.close()

    # (2) get the sequence lengths just for those proteins that we have domain information for:
    seq_lens, current_seq = {}, ''
    fasta_handle = gzip.open(seq_file, 'rt') if seq_file.endswith('gz') else open(seq_file)
    for fasta_line in fasta_handle:
        if fasta_line.startswith('>'):
            current_seq = fasta_line[1:-1].split()[0]
            seq_lens[current_seq] = 0
        else:
            seq_lens[current_seq] += len(fasta_line.strip())
    fasta_handle.close()

    # (3) write out the output!
    out_handle = gzip.open(out_trackfile, 'wt') if out_trackfile.endswith('gz') else open(out_trackfile, 'w')
    out_handle.write('# Pfam domain location tracks\n' +
                     '# Domains found in sequences from: '+seq_file+'\n' +
                     '# Original, complete domain hits found in: '+domain_locs+'\n' +
                     '\t'.join(['#Ensembl_Protein_ID', 'Domain_Name', 'Sequence_Length', 'Domain_Range'])+'\n')
    for seq_id in sorted(seq_lens.keys()):
        for domain, domain_range in sorted(list(prot_to_domains.get(seq_id, []))):
            out_handle.write(seq_id+'\t'+domain+'\t0-'+str(seq_lens[seq_id]-1)+'\t'+domain_range+'\n')
    out_handle.close()

    return out_trackfile


####################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process input PertInInt track files from Pfam domain hits.')
    parser.add_argument('--fasta_file', type=str,
                        help='Full path to a fasta-formatted file containing human protein sequences (for lengths)',
                        default=data_path + 'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.pep.all.fa')
    parser.add_argument('--domain_locations', type=str,
                        help='Full path to a tab-delimited file containing the domain locations in human proteins',
                        default=data_path + 'interacdome/Homo_sapiens.GRCh38.pep.all-domains-pfam_v31.tsv.gz')
    parser.add_argument('--outdir', type=str, default=data_path+'interacdome/',
                        help='Full path to write output track files to.')
    args = parser.parse_args()

    # create InteracDome input/output directory as needed:
    if not os.path.isdir(data_path+'interacdome'):
        call(['mkdir', data_path+'interacdome'])

    # ----------------------------------------------------------------------------------------------------
    if not os.path.isfile(args.fasta_file):
        sys.stderr.write('Could not find FASTA file in ' + args.fasta_file + '\n')
        sys.exit(1)

    # ----------------------------------------------------------------------------------------------------
    if not os.path.isfile(args.domain_locations):
        sys.stderr.write(
            'Could not find domain locations file in ' + args.domain_locations + '\n' +
            'Go to https://github.com/singh-lab/run-hmmer for complete instructions on creating this file.\n' +
            'Your input FASTA file should be: \n' +
            data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.pep.all.withgenelocs.verified.fa.gz\n'
        )
        sys.exit(1)

    # ----------------------------------------------------------------------------------------------------
    track_file = args.outdir + 'pfam'+PFAM_VERSION+'-domains_domsbyprot-' + GENOME_BUILD + '.txt.gz'
    if create_domain_tracks(args.fasta_file, args.domain_locations, track_file):
        sys.stderr.write('Wrote to '+track_file+'\n')
