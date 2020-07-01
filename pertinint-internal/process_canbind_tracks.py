#!/usr/bin/python

"""
Create track files corresponding to homologous regions to protein structural interfaces as determined
using some steps from Dario Ghersi's CanBind pipeline

Contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import gzip
import argparse
from config import data_path
from subprocess import call

####################################################################################################
# CONSTANTS
####################################################################################################

GENOME_BUILD = 'GRCh38'


####################################################################################################
# OBTAIN BIOLIP DATA
####################################################################################################

def download_initial_biolip_data(output_dir=data_path+'biolip/'):
    """
    :return: None, but download the initial release of BioLiP data (version 2013-03-6) if it has not yet
             been downloaded
    Download the most recent version of the BioLiP database. If you use this data, please cite:

    Jianyi Yang, Ambrish Roy and Yang Zhang (2013). "BioLiP: a semi-manually curated database for biologically relevant
    ligand-protein interactions." Nucleic Acids Research, 41: D1096-D1103.
    """

    for (bd_loc, basedata) in [
      (output_dir + 'downloaded_data/annotations/BioLiP.tar.bz2', 'BioLiP.tar.bz2'),
      (output_dir + 'downloaded_data/receptor1.tar.bz2', 'receptor1.tar.bz2'),
      (output_dir + 'downloaded_data/receptor_2013-03-13.tar.bz2', 'receptor_2013-03-13.tar.bz2'),
      (output_dir + 'downloaded_data/ligand_2013-03-6.tar.bz2', 'ligand_2013-03-6.tar.bz2')]:
        if not os.path.isfile(bd_loc):
            call(['wget', 'http://zhanglab.ccmb.med.umich.edu/BioLiP/download/' + basedata,
                  '-O', bd_loc])
        else:
            sys.stderr.write(bd_loc + ' has already been downloaded. Continuing...\n')

        # Expand the downloaded tar.bz2 files into the appropriate directories
        sys.stderr.write('Expanding ' + bd_loc + '\n')
        call(['tar', '-jxvf', bd_loc, '-C', '/'.join(bd_loc.split('/')[:-1]) + '/'])

    # note: for some reason, the receptor1 directory is missing several files from 2013-03-13.
    # we copy as many structures from receptor into receptor1 as possible to rescue some of these files
    for struct_file in os.listdir(output_dir + 'downloaded_data/receptor/'):
        start_with_one = False
        with open(output_dir+'downloaded_data/receptor/'+struct_file) as struct_handle:
            if struct_handle.next().strip().split()[5] == '1':
                start_with_one = True
        if start_with_one:
            sys.stderr.write('Copying '+struct_file+' to receptor1/ ...\n')
            call(['cp', output_dir + 'downloaded_data/receptor/'+struct_file,
                  output_dir + 'downloaded_data/receptor1/'+struct_file])

    # return the date that we have:
    return ['2013-03-6']


####################################################################################################

def update_biolip_data(output_dir=data_path+'biolip/'):
    """
    Download all available weekly updates to the original BioLiP data release
    :return: ordered list strings corresponding to all available weekly updates
    """

    weekly_updates_filename = output_dir + 'downloaded_data/weekly_updates.html'
    call(['wget', 'http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly.html', '-O', weekly_updates_filename])

    weekly_updates = []
    weekly_updates_handle = open(weekly_updates_filename)
    for html_line in weekly_updates_handle:
        # check if the date is contained on this line
        if html_line.startswith('<tr><td>') and html_line.strip().endswith('</td>') and \
           len(html_line.strip()) in [22, 23]:
            weekly_updates.append(html_line.strip().replace('<tr><td>', '').replace('</td>', ''))
    weekly_updates_handle.close()
    call(['rm', weekly_updates_filename])

    # Now, download and expand each weekly update that we haven't already downloaded
    for date in weekly_updates:
        if not os.path.isfile(output_dir + 'downloaded_data/annotations/BioLiP_' + date + '.txt'):
            call(['wget', 'http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly/BioLiP_' + date + '.txt',
                  '-O', output_dir + 'downloaded_data/annotations/BioLiP_' + date + '.txt'])

            for tarball_type in ['ligand', 'receptor1']:
                call(['wget',
                      'http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly/' + tarball_type + '_' + date + '.tar.bz2',
                      '-O', output_dir + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2'])
                call(['tar', '-jxvf', output_dir + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2',
                      '-C', output_dir + 'downloaded_data/'])
                call(['rm', output_dir + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2'])

    sys.stderr.write('All BioLiP data through ' + weekly_updates[0] + ' has been downloaded!\n')

    return weekly_updates + ['2013-03-6']  # include the base annotation date


####################################################################################################

def create_biolip_annotation_file(output_dir, release_dates):
    """
    :param output_dir: full path to a directory where BioLiP data has been downloaded
    :param release_dates: ordered list of "release dates" from BioLiP for which we have downloaded data
    :return: none, but print error message if we could not create an up-to-date annotation file
    """

    # No matter what, create a current concatenated list of *all* annotations
    current_annotations = [output_dir + 'downloaded_data/annotations/BioLiP_' + update + '.txt'
                           for update in release_dates
                           if os.path.isfile(output_dir + 'downloaded_data/annotations/BioLiP_' + update + '.txt')]
    failed_annotations = [output_dir + 'downloaded_data/annotations/BioLiP_' + update + '.txt'
                          for update in release_dates
                          if not os.path.isfile(output_dir + 'downloaded_data/annotations/BioLiP_' + update + '.txt')]

    concatenated_annotation_file = output_dir + 'processed_data/annotations/current_annotations.txt'
    file_handle = open(concatenated_annotation_file, 'w')
    call(['cat'] + current_annotations, stdout=file_handle)
    file_handle.close()

    if len(failed_annotations) > 0:
        sys.stderr.write('Failed to update results from the following BioLiP releases:\n' +
                         '\n'.join(failed_annotations) + '\n')


####################################################################################################

def flag_missing_structures(biolip_dir=data_path+'biolip/'):
    """
    :param annotation_file: full path to a concatenated annotation file
    :param annotation_path: full path to a list of annotations by weekly update
    :param ligand_path: full path to a directory containing all ligand structures
    :param receptor_path: full path to a directory containing all receptor structures (renumbered to start at 1)
    :return: list of receptor and ligand structures that should be present but are not
    """

    annotation_file = biolip_dir + 'processed_data/annotations/current_annotations.txt'
    ligand_pdb_dir = biolip_dir + 'downloaded_data/ligand/'
    receptor_pdb_dir = biolip_dir + 'downloaded_data/receptor1/'

    failure_file = biolip_dir + 'processed_data/annotations/failed_ligand_annotations.txt'
    failure_ligand_handle = gzip.open(failure_file, 'wt') if failure_file.endswith('gz') else open(failure_file, 'w')
    failure_file = biolip_dir + 'processed_data/annotations/failed_receptor_annotations.txt'
    failure_receptor_handle = gzip.open(failure_file, 'wt') if failure_file.endswith('gz') else open(failure_file, 'w')

    # go through the annotation file and keep track of structures that are missing:
    missing_structures = set()
    annot_handle = gzip.open(annotation_file) if annotation_file.endswith('gz') else open(annotation_file)
    for annot_line in annot_handle:
        (annot_pdb_id, pdb_chain, _, _, ligand_id, ligand_chain, ligand_serial_number) = annot_line[:-1].split('\t')[:7]

        ligand_structure_file = ligand_pdb_dir + annot_pdb_id + '_' + ligand_id + '_' + ligand_chain + '_' + \
                                ligand_serial_number + '.pdb'
        receptor_structure_file = receptor_pdb_dir + annot_pdb_id + pdb_chain + '.pdb'

        if not os.path.isfile(ligand_structure_file):
            failure_ligand_handle.write(annot_line)
            missing_structures.add(annot_pdb_id)
        if not os.path.isfile(receptor_structure_file):
            failure_receptor_handle.write(annot_line)
            missing_structures.add(annot_pdb_id)
    annot_handle.close()
    failure_ligand_handle.close()
    failure_receptor_handle.close()

    # determine which dates need to be redownloaded
    if len(missing_structures) > 0:
        annotation_dir = biolip_dir + 'downloaded_data/annotations/'
        missing_dates = set()

        for annot_file in [annotation_dir+f for f in os.listdir(annotation_dir) if f.endswith('.txt')]:
            date = annot_file.split('/')[-1][7:-4]  # e.g., "2013-03-13" from "BioLiP_2013-03-13.txt"
            with open(annot_file) as annot_handle:
                for annot_line in annot_handle:
                    if annot_line[:-1].split('\t')[0] in missing_structures:
                        missing_dates.add(date)
                        break

        # redownload those dates that are missing!
        download_new_data = False
        if len(missing_dates) > 0:
            for date in sorted(list(missing_dates)):
                sys.stderr.write(date+'\n')
                if download_new_data:
                    call(['wget', 'http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly/BioLiP_' + date + '.txt',
                          '-O', biolip_dir + 'downloaded_data/annotations/BioLiP_' + date + '.txt'])

                    for tarball_type in ['ligand', 'receptor1']:
                        call(['wget',
                              'http://zhanglab.ccmb.med.umich.edu/BioLiP/weekly/' + tarball_type + '_' + date + '.tar.bz2',
                              '-O', biolip_dir + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2'])
                        call(['tar', '-jxvf',
                              biolip_dir + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2',
                              '-C', biolip_dir + 'downloaded_data/'])
                        call(['rm', biolip_dir + 'downloaded_data/' + tarball_type + '_' + date + '.tar.bz2'])


####################################################################################################
# SPLIT CLUSTERS BY BINDING TYPE
####################################################################################################

def find_ions(lig_grp=data_path+'interacdome/downloaded_data/ligand_groups.txt'):
    """
    :param lig_grp:
    :return:
    """

    if not os.path.isfile(lig_grp):
        call(['wget',
              'https://raw.githubusercontent.com/Singh-Lab/InteracDome/master/downloaded_data/ligand_groups.txt',
              '-O', lig_grp])

    ions = set()
    lig_handle = gzip.open(lig_grp, 'rt') if lig_grp.endswith('gz') else open(lig_grp)
    for lig_line in lig_handle:
        if lig_line.startswith('ION_'):
            ions.add(lig_line.strip().split()[-1])

    return ions


####################################################################################################

def reformat_cluster_files(ion_set, binding_results_dir=data_path+'canbind/pipeline/step.03/' +
                                                        'bindingInfoBioLip_ID60_COV80_BINDID90/'):
    """
    :param ion_set: set of ligands that would be considered ions
    :param binding_results_dir: full path to a directory containing clusters and structures as computed
                                using CanBind.
    :return: none
    """

    for cluster_file in [binding_results_dir+fname for fname in os.listdir(binding_results_dir)
                         if fname.endswith('.clusters')]:
        clust_handle = open(cluster_file)
        cluster_out = open(cluster_file+'-new', 'w')

        for clustered_bindsites in clust_handle:  # each line is a new binding site!

            # separate each match into one of the following categories:
            matches_by_ligtype = {}

            for match in clustered_bindsites[:-1].split(' - '):
                c = match.strip().split('\t')
                if len(c) < 5:
                    continue
                structure_id, orig_bindsite_name, orig_bindpos, updated_bindpos, ligand_type = c[:5]

                # get the ligands
                if ligand_type == 'III':
                    grp = 'III'
                elif ligand_type == 'NUC':
                    grp = 'NUC'
                elif ligand_type in ion_set:
                    grp = 'ION'
                else:
                    grp = 'SM'

                if grp not in matches_by_ligtype:
                    matches_by_ligtype[grp] = []
                matches_by_ligtype[grp].append(match)

            for matches in matches_by_ligtype.values():
                cluster_out.write(' - '.join(matches)+'\n')
        clust_handle.close()
        cluster_out.close()

        call(['mv', cluster_file+'-new', cluster_file])


####################################################################################################
# CREATE TRACKS
####################################################################################################

def subtype_nuc_ligand(ligand_structure_file):
    """
    :param ligand_structure_file: full path to a ligand structural PDB file (downloaded from BioLiP)
    :return: set of ligands (RNA, DNA) found in the file
    """

    # check that the file exists and that it is of type "nucleic acid"
    if not os.path.isfile(ligand_structure_file) or '_NUC_' not in ligand_structure_file:
        return set()

    # keep track of ligand types observed:
    # ligand_types = set()

    curr_base_id = ''  # keep track of the current base ID
    curr_base_lines = []  # and all lines from the original ligand file corresponding to that base

    orig_ligand_handle = open(ligand_structure_file)
    for lig_line in orig_ligand_handle:

        # inappropriately short / non-ATOM entry line:
        if len(lig_line) < 81:
            continue

        # reached a new base, so process and reset:
        if curr_base_id != lig_line[22:26].strip():

            if len(curr_base_lines) > 0:

                # "deoxy-ribonucleic acid" (DNA) will not have any O2 backbone atoms
                atoms = [base_line[12:16].strip() for base_line in curr_base_lines]
                if "O2'" in atoms:
                    return 'RNA'
                # ligand_types.add('RNA' if "O2'" in atoms else 'DNA')

            # reset the current list of lines for the bext base to process:
            curr_base_lines = []
            curr_base_id = lig_line[22:26].strip()

        # no matter what, add the current line (to be eventually processed)
        curr_base_lines.append(lig_line)
    orig_ligand_handle.close()

    # repeat the processing step for the very last set of entries:
    if len(curr_base_lines) > 0:
        atoms = [base_line[12:16].strip() for base_line in curr_base_lines]
        if "O2'" in atoms:
            return 'RNA'
        # ligand_types.add('RNA' if "O2'" in atoms else 'DNA')

    return 'DNA'
    # return ligand_types  # return the set of ligand types found here


####################################################################################################

def create_rna_structure_list(biolip_dir, rna_output):
    """
    :param biolip_dir: full path to the biolip directory containing a compiled list of current annotations as
                       well as a directory with all ligand PDBs
    :param rna_output: full path to a file to write a list of PDB identifiers out to if they were found in complex
                       with RNA
    :return: none
    """

    rna_structures = set()  # set of PDB IDs that are found in complex with RNA

    annotation_file = biolip_dir + 'processed_data/annotations/current_annotations.txt'
    ligand_pdb_dir = biolip_dir + 'downloaded_data/ligand/'

    # go through the annotation file and reanalyze ligand PDB files labeled as "NUC"
    annot_handle = gzip.open(annotation_file) if annotation_file.endswith('gz') else open(annotation_file)
    for annot_line in annot_handle:
        (annot_pdb_id, pdb_chain, _, _, ligand_id, ligand_chain, ligand_serial_number) = annot_line[:-1].split('\t')[:7]

        if ligand_id == 'NUC' and annot_pdb_id not in rna_structures:  # if this PDB ID is already found w/ RNA, skip
            ligand_structure_file = ligand_pdb_dir + annot_pdb_id + '_' + ligand_id + '_' + ligand_chain + '_' + \
                                    ligand_serial_number + '.pdb'
            ligand_types = subtype_nuc_ligand(ligand_structure_file)
            if ligand_types == 'RNA':  # if 'RNA' in ligand_types:
                rna_structures.add(annot_pdb_id)
    annot_handle.close()

    # write out to the RNA file all PDB IDs that correspond to RNA (rather than or in addition to DNA)
    rna_handle = gzip.open(rna_output, 'wt') if rna_output.endswith('gz') else open(rna_output, 'w')
    rna_handle.write('\n'.join(sorted(list(rna_structures)))+'\n')
    rna_handle.close()


####################################################################################################

def create_homology_tracks(out_trackfile, out_weightfile, seq_file, biolip_dir,
                           rna_structures_file=data_path+'canbind/pipeline/step.03/rna_structures.txt',
                           binding_results_dir=data_path+'canbind/pipeline/step.03/' +
                                               'bindingInfoBioLip_ID60_COV80_BINDID90/'):
    """
    :param out_trackfile: full path to a file to write out track information
    :param out_weightfile: full path to a file to write out per-position track weights
    :param seq_file: full path to the FASTA-formatted sequence file used by CanBind
    :param biolip_file: full path to the tab-delimited BioLiP annotation file reflecting BioLiP's update date
    :param rna_structures_file: full path to a file containing a list of PDB structures that contain RNA
                                as we currently cannot distinguish DNA from RNA
    :param binding_results_dir: full path to a directory containing clusters and structures as computed
                                using CanBind.
    :return: None, but print success message upon write of two output files
    """

    # (1) Get the last date updated for the BioLiP file and the number of sequences in the FASTA file:
    annot_dates = set()
    for annot_file in os.listdir(biolip_dir+'downloaded_data/annotations/'):
        if annot_file.startswith('BioLiP_') and annot_file.endswith('.txt') and '-' in annot_file:
            annot_dates.add(tuple(map(int, annot_file[7:-4].split('-'))))
    biolip_update_date = '-'.join([str(date) for date in sorted(list(annot_dates), reverse=True)[0]])

    seq_count = 0
    fasta_handle = gzip.open(seq_file) if seq_file.endswith('gz') else open(seq_file)
    for fasta_line in fasta_handle:
        if fasta_line.startswith('>'):
            seq_count += 1
    fasta_handle.close()

    # (2) write header to track file
    domlocs_handle = gzip.open(out_trackfile, 'wt') if out_trackfile.endswith('gz') else open(out_trackfile, 'w')
    domlocs_handle.write(
        '# All distinct binding sites that were found by mapping structural information onto sequences by CanBind\n')
    domlocs_handle.write(
        '# CanBind was trained on all BioLiP sequences and structures available as of '+biolip_update_date+', ' +
        'and was run on all '+"{:,}".format(seq_count)+' sequences in\n# '+seq_file+'\n')
    domlocs_handle.write('# Original weights can be found in ' + binding_results_dir + '\n')
    domlocs_handle.write('\t'.join(['#EnsemblProtID', 'CanBindBindingSite', 'matchstate:AA-0-index:AA-value']) + '\n')

    # (3) write header to weights file
    weights_handle = gzip.open(out_weightfile, 'wt') if out_weightfile.endswith('gz') else open(out_weightfile, 'w')
    weights_handle.write('# Continuous weights across inferred binding site residues from CanBind, by ligand type\n')
    weights_handle.write(
        '# CanBind was trained on all BioLiP sequences and structures available as of '+biolip_update_date+', ' +
        'and was run on all '+"{:,}".format(seq_count)+' sequences in\n# '+seq_file+'\n')
    weights_handle.write(
        '\t'.join(['#EnsemblProtID_BindingSiteNum', 'LigandType', '1-Indexed-MatchState', 'FracWithin4',
                   'NumInstancesAlways10', 'FullMatchState']) + '\n')

    # (4) Get the list of RNA structures (to separate accordingly)
    rna_handle = gzip.open(rna_structures_file) if rna_structures_file.endswith('gz') else open(rna_structures_file)
    rna_structures = set([rna_line.strip() for rna_line in rna_handle])
    rna_handle.close()

    # (5) get the list of proteins that we have some sort of homology mapping of interactions to:
    prots = sorted(list(set([a.replace('.binding', '').replace('.clusters', '').replace('.struct', '')
                             for a in os.listdir(binding_results_dir) if not os.path.isdir(binding_results_dir + a)])))

    # (6) process each protein:
    for prot_id in prots:
        if os.path.isfile(binding_results_dir + prot_id + '.clusters'):
            binding_file = binding_results_dir + prot_id + '.clusters'
        else:
            binding_file = binding_results_dir + prot_id + '.binding'

        with open(binding_file) as binding_handle:
            binding_sites = [bind_line.strip() for bind_line in binding_handle]
        with open(binding_results_dir + 'weights/' + prot_id + '.weighted') as fracweights_handle:
            all_weights = [fracweights_line[:-1] for fracweights_line in fracweights_handle]

        for bindsite_id, clustered_bindsites in enumerate(binding_sites):  # each line is a new binding site!

            combined_matches = [c.strip().split('\t') for c in clustered_bindsites.split('-')]
            current_ligands = {}  # keep track of all ligands seen (& their COUNTS to order them in case of clustering)
            allindices = set()  # keep track of the 0-indices of this binding site, too

            for c in combined_matches:
                if len(c) < 5:
                    continue
                structure_id, orig_bindsite_name, orig_bindpos, updated_bindpos, ligand_type = c[:5]

                # get the ligands
                cltype = ligand_type
                if ligand_type == 'III':
                    cltype = 'PEPTIDE_'
                elif ligand_type == 'NUC':
                    cltype = 'RNA_' if structure_id[:4] in rna_structures else 'DNA_'

                if cltype not in current_ligands:
                    current_ligands[cltype] = 0
                current_ligands[cltype] += 1

                # get the binding sites
                for m in updated_bindpos.split():
                    allindices.add(int(m[1:]) - 1)

            # list of all ligands:
            all_ligands = ','.join([lig for ligcnt, lig in sorted([(v, k) for k, v in current_ligands.items()],
                                                                  reverse=True)])

            # Now, let's get the complete "domain" (we'll worry about the weights later!)
            thisseq = [seqval.split('\t')[1] for seqval in open(binding_results_dir + prot_id + '.struct')]
            thisdom = range(min(allindices), max(allindices) + 1)

            mstates = []
            for matchstate, dompos in enumerate(thisdom):
                mstates.append(str(matchstate + 1) + ':' + str(dompos) + ':' + thisseq[dompos])
            domlocs_handle.write(
                '\t'.join([prot_id, prot_id + '_BS' + str(bindsite_id + 1).zfill(2), ','.join(mstates)]) + '\n')

            # AND, let's write out the appropriate weights...
            nonzeros = all_weights[bindsite_id].split('\t')
            for m in nonzeros:
                dompos = int(m.split()[0][1:]) - 1
                domval = m.split()[0][0]
                weight = m.split()[1]

                weights_handle.write('\t'.join([prot_id + '_BS' + str(bindsite_id + 1).zfill(2), all_ligands,
                                                str(thisdom.index(dompos) + 1),
                                                weight, '10',
                                                str(thisdom.index(dompos) + 1) + ':' + domval]) + '\n')
    weights_handle.close()
    domlocs_handle.close()


####################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Process input PertInInt track files from homology-based binding ' +
                                                 'regions obtained using CanBind.')
    parser.add_argument('--initialize_biolip', dest='initialize_biolip', action='store_true', default=False,
                        help='Download the starting set of structures and annotations (from March 6, 2013)')
    parser.add_argument('--update_biolip', dest='update_biolip', action='store_true', default=False,
                        help='Update the downloaded structures from BioLiP with the most current version.')
    parser.add_argument('--update_failed_biolip', dest='update_failed_biolip', action='store_true', default=False,
                        help='Determine which updates failed to download and extract properly; redownload as needed.')
    parser.add_argument('--fix_clusters', dest='fix_clusters', action='store_true', default=False,
                        help='Split clusters formed by CanBind such that no single cluster contains ligands from ' +
                             'multiple ligand "classes" (i.e., nucleic acids, peptides, ions, small molecules)')
    parser.add_argument('--biolip_dir', type=str, default=data_path+'biolip/',
                        help='Full path to directory containing structures downloaded from BioLiP.')
    parser.add_argument('--fasta_file', type=str,
                        default=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.pep.all.fa',
                        help='Full path to the FASTA file of human sequences against which we searched for ' +
                             'homologous interaction sites.')
    parser.add_argument('--rna_file', type=str, default=data_path+'canbind/pipeline/data/biolip/rna_structure_list.txt',
                        help='Full path to a list of PDB IDs that are found in complex with RNA.')
    parser.add_argument('--outdir', type=str, default=data_path+'canbind/',
                        help='Full path to write output track files to.')
    args = parser.parse_args()

    # ----------------------------------------------------------------------------------------------------
    # Download the initial set of structures and annotations, *if need be*
    if args.initialize_biolip:
        sys.stderr.write('Downloading original BioLiP data...\n')

        # create CanBind input/output directory as needed:
        for subdirectory in [args.biolip_dir,
                             args.biolip_dir + 'downloaded_data',
                             args.biolip_dir + 'downloaded_data/annotations',
                             args.biolip_dir + 'downloaded_data/ligand',
                             args.biolip_dir + 'downloaded_data/receptor1',
                             args.biolip_dir + 'processed_data',
                             args.biolip_dir + 'processed_data/annotations']:
            if not os.path.isdir(subdirectory):
                call(['mkdir', subdirectory])

        # Check to see if starting tar.bz2 files already exist; download and expand otherwise
        release_dates = download_initial_biolip_data(args.biolip_dir)

        # Create a current concatenated list of *all* annotations
        create_biolip_annotation_file(args.biolip_dir, release_dates)

    # ----------------------------------------------------------------------------------------------------
    elif args.update_biolip:
        sys.stderr.write('Updating BioLiP data...\n')

        # Make sure that the initialization step has already been run
        for subdirectory in [args.biolip_dir,
                             args.biolip_dir + 'downloaded_data',
                             args.biolip_dir + 'downloaded_data/annotations',
                             args.biolip_dir + 'downloaded_data/ligand',
                             args.biolip_dir + 'downloaded_data/receptor1',
                             args.biolip_dir + 'processed_data',
                             args.biolip_dir + 'processed_data/annotations']:
            if not os.path.isdir(subdirectory):
                sys.stderr.write('No such directory: ' + subdirectory + '\n')
                sys.stderr.write('Please run python ' + sys.argv[0] + ' --initialize_biolip ' +
                                 '--biolip_dir '+args.biolip_dir+'\n')
                sys.exit(1)

        # Get the list of all weekly updates (this does not include the base data) and download corresponding files
        release_dates = update_biolip_data(args.biolip_dir)

        # Create a current concatenated list of *all* annotations
        create_biolip_annotation_file(args.biolip_dir, release_dates)

    # ----------------------------------------------------------------------------------------------------
    elif args.update_failed_biolip:
        sys.stderr.write('Checking for missing BioLiP downloads...\n')
        flag_missing_structures(args.biolip_dir)

    # ----------------------------------------------------------------------------------------------------
    elif args.fix_clusters:
        sys.stderr.write('Confirming that clusters do not contain multiple ligand classes...\n')
        ion_set = find_ions(data_path+'canbind/pipeline/step.03/ligand_groups.txt')
        reformat_cluster_files(ion_set)

    else:
        if not os.path.isfile(args.fasta_file):
            sys.stderr.write('Could not find FASTA file in ' + args.fasta_file + '\n')
            sys.exit(1)

        # ------------------------------------------------------------------------------------------------
        if not os.path.isdir(args.biolip_dir):
            sys.stderr.write(
                'Could not find BioLiP downloads in ' + args.biolip_dir + '. Please run:\n' +
                'python ' + sys.argv[0] + ' --initialize_biolip --biolip_dir '+args.biolip_dir+'\n',
                'python ' + sys.argv[0] + ' --update_biolip --biolip_dir ' + args.biolip_dir + '\n'
            )
            sys.exit(1)

        # ------------------------------------------------------------------------------------------------
        if not os.path.isfile(args.rna_file):
            sys.stderr.write('Finding ligands that contain RNA... ')
            create_rna_structure_list(args.biolip_dir, args.rna_file)
            sys.stderr.write('Done!\n')

        # ------------------------------------------------------------------------------------------------
        track_file = args.outdir + 'canbind-biolip-to-ensembl_domsbyprot-' + GENOME_BUILD + '.txt.gz'
        weight_file = args.outdir + 'canbind-biolip-to-ensembl_domainweights-' + GENOME_BUILD + '.txt.gz'
        create_homology_tracks(track_file, weight_file, args.fasta_file, args.biolip_dir, args.rna_file)
