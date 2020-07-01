#!/usr/bin/python

"""
Create necessary track files reflecting interaction-based functionality weights from InteracDome

Contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import argparse
import gzip
from config import data_path, GENOME_BUILD
from subprocess import call

####################################################################################################
# CONSTANTS
####################################################################################################

PFAM_VERSION = '31'


####################################################################################################

def get_restricted_domains(minimum_instances, minimum_structures, stricter_instance_cutoff,
                           id_file=data_path+'interacdome/InteracDome_v0.3-confident.tsv'):
    """
    :param minimum_instances: minimum number of instances that a (domain, ligand) pair must have to be considered
    :param minimum_structures: minimum number of structures that a (domain, ligand) pair must have to be considered
    :param id_file: full path to a list of confidently modeled interaction domains downloaded from InteracDome
    :return: none, but reset RESTRICTED_DOMS accordingly...
    """

    skip_doms = set()

    # remove the common false positive repetitive domain/small molecule groupings:
    for problem_domain in ['PF07686_V-set', 'PF00047_ig', 'PF07679_I-set', 'PF13927_Ig_3', 'PF13895_Ig_2',
                           'PF08205_C2-set_2', 'PF07654_C1-set', 'PF00008_EGF']:
        for problem_ligand in ['BMA', 'MAN', 'UNL', 'GAL', 'FUC']:
            skip_doms.add((problem_domain, problem_ligand))

    # remove DNA "mimicry" by RNA structures:
    for rna_problem_domain in ['PF00047_ig', 'PF00096_zf-C2H2']:
        for rna in ['RNA_', 'RNABASE_', 'RNABACKBONE_']:
            skip_doms.add((rna_problem_domain, rna))

    # and remove RNA "mimicry" by DNA structures:
    for dna_problem_domain in ['PF04851_ResIII']:
        for dna in ['DNA_', 'DNABASE_', 'DNABACKBONE_']:
            skip_doms.add((dna_problem_domain, dna))

    # finally, remove domains that do not have the required number of instances/structures in the BioLiP
    with open(id_file) as id_handle:
        header = None
        for idline in id_handle:
            if idline.startswith('#'):
                continue
            elif not header:
                header = idline[:-1].split('\t')
                continue
            v = idline[:-1].split('\t')

            # any ligand type that doesn't pass the bare minimum threshold will be skipped
            if int(v[header.index('num_nonidentical_instances')]) < minimum_instances or \
               int(v[header.index('num_structures')]) < minimum_structures:
                skip_doms.add((v[header.index('pfam_id')], v[header.index('ligand_type')]))

            # we require an even higher threshold for small molcules and specific ions:
            if v[header.index('ligand_type')] not in {'DNA_', 'DNABASE_', 'DNABACKBONE',
                                                      'RNA_', 'RNABASE_', 'RNABACKBONE_',
                                                      'PEPTIDE_'} and \
               int(v[header.index('num_nonidentical_instances')]) < stricter_instance_cutoff:
                skip_doms.add((v[header.index('pfam_id')], v[header.index('ligand_type')]))

    return skip_doms


####################################################################################################

def create_interacdome_tracks(interacdome_wts, domain_locs, skip_doms, out_trackfile, out_weightfile):
    """
    :param interacdome_wts: full path to a tab-delimited file downloaded from InteracDome containing
                               per-position domain binding frequencies by ligand type
    :param domain_locs: full path to a tab-delimited file containing domain locations in human proteins
    :param skip_doms: set of (domain name, ligand type) pairs to *not* include
    :param out_trackfile: full path to a file to write out track information
    :param out_weightfile: full path to a file to write out per-position track weights
    :return: None, but print success message upon write of two output files
    """

    interaction_domains = set()  # keep track of Pfam domain identifiers that we have binding information for

    # (1) reformat the InteracDome results to work with PertInInt
    weight_handle = gzip.open(out_weightfile, 'wt') if out_weightfile.endswith('gz') else open(out_weightfile, 'w')

    interacdome_handle = gzip.open(interacdome_wts, 'rt') if interacdome_wts.endswith('gz') else open(interacdome_wts)
    header = None
    for wt_line in interacdome_handle:
        if wt_line.startswith('#'):
            weight_handle.write(wt_line)
            continue
        elif not header:
            header = wt_line[:-1].split('\t')
            weight_handle.write('# Scores REFORMATTED for use by PertInInt\n' +
                                '\t'.join(['#domain_name', 'ligand_type', '1-indexed_match_state',
                                           'binding_frequency', 'confidence_always10', 'misc'])+'\n')
            continue

        values = wt_line[:-1].split('\t')
        domain_name = values[header.index('pfam_id')]
        ligand_type = values[header.index('ligand_type')]

        if (domain_name, ligand_type) in skip_doms:
            continue

        binding_freqs = [float(v) for v in values[header.index('binding_frequencies')].split(',')]
        for matchstate, bf in enumerate(binding_freqs):
            if bf > 0:
                interaction_domains.add(domain_name)
                weight_handle.write('\t'.join([domain_name, ligand_type, str(matchstate+1), str(bf), '10',
                                              str(matchstate+1)+':X'])+'\n')
    interacdome_handle.close()
    weight_handle.close()

    # (2) turns out the domain file is already in the format we want! we'll restrict to those domains
    #     that we have interaction data for
    track_handle = gzip.open(out_trackfile, 'wt') if out_trackfile.endswith('gz') else open(out_trackfile, 'w')
    track_handle.write('# Domain hits LIMITED to those with InteracDome scores ('+interacdome_wts+')\n')

    domloc_handle = gzip.open(domain_locs, 'rt') if domain_locs.endswith('gz') else open(domain_locs)
    for dom_line in domloc_handle:
        if dom_line.startswith('#'):
            track_handle.write(dom_line)
            continue

        # write those lines corresponding to interaction domains
        pfam_id = dom_line[:-1].split('\t')[1]
        if pfam_id in interaction_domains:
            track_handle.write(dom_line)
    domloc_handle.close()
    track_handle.close()


####################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process input PertInInt track files from InteracDome scores.')
    parser.add_argument('--interacdome_scores', type=str,
                        help='Full path to a tab-delimited file containing the InteracDome binding values.',
                        default=data_path + 'interacdome/InteracDome_v0.3-confident.tsv')
    parser.add_argument('--domain_locations', type=str,
                        help='Full path to a tab-delimited file containing the domain locations in human proteins',
                        default=data_path + 'interacdome/Homo_sapiens.GRCh38.pep.all.withgenelocs.verified' +
                                '-domains-pfam_v31.tsv.gz')
    parser.add_argument('--outdir', type=str, default=data_path+'interacdome/',
                        help='Full path to write output track files to.')
    args = parser.parse_args()

    # create InteracDome input/output directory as needed:
    if not os.path.isdir(data_path+'interacdome'):
        call(['mkdir', data_path+'interacdome'])

    # ----------------------------------------------------------------------------------------------------
    if not os.path.isfile(args.interacdome_scores):
        sys.stderr.write(
            'Could not find InteracDome file in ' + args.interacdome_scores + '\n' +
            'Go to https://interacdome.princeton.edu > Download > Confident Domain-Ligand Interactions.\n' +
            'Save this file to '+args.interacdome_scores+'\n'
        )
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
    restricted_doms = get_restricted_domains(5., 0., 7.5, args.interacdome_scores)
    track_file = args.outdir + 'interacdome0.3-pfam'+PFAM_VERSION+'_domsbyprot-' + GENOME_BUILD + '.txt.gz'
    weight_file = args.outdir + 'interacdome0.3-pfam'+PFAM_VERSION+'_domainweights-' + GENOME_BUILD + '.txt.gz'
    create_interacdome_tracks(args.interacdome_scores, args.domain_locations, restricted_doms, track_file, weight_file)
