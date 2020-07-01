#!/usr/bin/python

"""
Separate the 100-way vertebrate alignment(s) from UCSC's Genome Browser into separate directories
corresponding to the known protein-coding genes

Contact Shilpa N. Kobren (snadimpa@alumni.princeton.edu) with questions
"""

import os
import sys
import argparse
import gzip
import math
from subprocess import call
from Bio import AlignIO
from intervaltree import IntervalTree
from config import data_path, GENOME_BUILD, BUILD_ALT_ID

####################################################################################################
# CONSTANTS
####################################################################################################

if GENOME_BUILD == 'GRCh38':
    SPCS = ['GRCh38/hg38', 'CSAC 2.1.4/panTro4', 'gorGor3.1/gorGor3', 'WUGSC 2.0.2/ponAbe2', 'GGSC Nleu3.0/nomLeu3',
            'BGI CR_1.0/rheMac3', 'Macaca_fascicularis_5.0/macFas5', 'Baylor Panu_2.0/papAnu2',
            'Chlorocebus_sabeus 1.1/chlSab2', 'WUGSC 3.2/calJac3', 'Broad/saiBol1', 'Broad/otoGar3',
            'TupChi_1.0/tupChi1', 'Broad/speTri2', 'JacJac1.0/jacJac1', 'MicOch1.0/micOch1',
            'C_griseus_v1.0/criGri1', 'MesAur1.0/mesAur1', 'GRCm38/mm10', 'RGSC 6.0/rn6',
            'Broad HetGla_female_1.0/hetGla2', 'Broad/cavPor3', 'ChiLan1.0/chiLan1', 'OctDeg1.0/octDeg1',
            'Broad/oryCun2', 'OchPri3.0/ochPri3', 'SGSC Sscrofa10.2/susScr3', 'Vicugna_pacos-2.0.1/vicPac2',
            'CB1/camFer1',
            'Baylor Ttru_1.4/turTru2', 'Oorc_1.1/orcOrc1', 'PHO1.0/panHod1', 'UMD_3.1.1/bosTau8',
            'ISGC Oar_v3.1/oviAri3',
            'CHIR_1.0/capHir1', 'Broad/equCab2', 'CerSimSim1.0/cerSim1', 'ICGSC Felis_catus 8.0/felCat8',
            'Broad CanFam3.1/canFam3', 'MusPutFur1.0/musFur1', 'BGI-Shenzhen 1.0/ailMel1', 'Oros_1.0/odoRosDiv1',
            'LepWed1.0/lepWed1', 'ASM32557v1/pteAle1', 'Broad/pteVam1', 'ASM32734v1/myoDav1',
            'Broad Institute Myoluc2.0/myoLuc2', 'EptFus1.0/eptFus1', 'EriEur2.0/eriEur2', 'Broad/sorAra2',
            'ConCri1.0/conCri1', 'Broad/loxAfr3', 'EleEdw1.0/eleEdw1', 'Broad v1.0/triMan1', 'ChrAsi1.0/chrAsi1',
            'Broad/echTel2', 'OryAfe1.0/oryAfe1', 'Baylor/dasNov3', 'Broad/monDom5', 'WTSI Devil_ref v7.0/sarHar1',
            'TWGS Meug_1.1/macEug2', 'WUGSC 5.0.1/ornAna1', 'F_cherrug_v1.0/falChe1', 'F_peregrinus_v1.0/falPer1',
            'FicAlb1.5/ficAlb2', 'ASM38545v1/zonAlb1', 'GeoFor_1.0/geoFor1', 'WashU taeGut324/taeGut2',
            'PseHum1.0/pseHum1',
            'WUSTL v6.3/melUnd1', 'AV1/amaVit1', 'SMACv1.1/araMac1', 'Cliv_1.0/colLiv1', 'BGI_duck_1.0/anaPla1',
            'ICGSC Gallus_gallus-4.0/galGal4', 'TGC Turkey_2.01/melGal1', 'allMis0.2/allMis1', 'CheMyd_1.0/cheMyd1',
            'v3.0.3/chrPic2', 'PelSin_1.0/pelSin1', 'ASM38561v1/apaSpi1', 'Broad AnoCar2.0/anoCar2', 'JGI 7.0/xenTro7',
            'Broad/latCha1', 'Genoscope 8.0/tetNig2', 'FUGU5/fr3', 'version 1 of Takifugu flavidus genome/takFla1',
            'Broad oreNil1.1/oreNil2', 'NeoBri1.0/neoBri1', 'AstBur1.0/hapBur1', 'MetZeb1.1/mayZeb1',
            'PunNye1.0/punNye1',
            'NIG/UT MEDAKA1/oryLat2', 'Xiphophorus_maculatus-4.4.2/xipMac1', 'Broad/gasAcu1',
            'Genofisk GadMor_May2010/gadMor1', 'Zv10/danRer10', 'Astyanax_mexicanus-1.0.2/astMex1',
            'LepOcu1/lepOcu1', 'WUGSC 7.0/petMar2']
elif GENOME_BUILD == 'GRCh37':
    SPCS = ['GRCh37/hg19', 'CSAC 2.1.4/panTro4', 'gorGor3.1/gorGor3', 'WUGSC 2.0.2/ponAbe2', 'GGSC Nleu3.0/nomLeu3',
            'BGI CR_1.0/rheMac3', 'Macaca_fascicularis_5.0/macFas5', 'Baylor Pham_1.0/papHam1',
            'Chlorocebus_sabeus 1.0/chlSab1', 'WUGSC 3.2/calJac3', 'Broad/saiBol1', 'Broad/otoGar3',
            # Euarchontoglires subset
            'TupChi_1.0/tupChi1', 'Broad/speTri2', 'JacJac1.0/jacJac1', 'MicOch1.0/micOch1', 'C_griseus_v1.0/criGri1',
            'MesAur1.0/mesAur1', 'GRCm38/mm10', 'RGSC 5.0/rn5', 'Broad HetGla_female_1.0/hetGla2', 'Broad/cavPor3',
            'ChiLan1.0/chiLan1', 'OctDeg1.0/octDeg1', 'Broad/oryCun2', 'OchPri3.0/ochPri3',
            # Laurasiatheria subset
            'SGSC Sscrofa10.2/susScr3', 'Vicugna_pacos-2.0.1/vicPac2', 'CB1/camFer1', 'Baylor Ttru_1.4/turTru2',
            'Oorc_1.1/orcOrc1', 'PHO1.0/panHod1', 'Baylor Btau_4.6.1/bosTau7', 'ISGC Oar_v3.1/oviAri3',
            'CHIR_1.0/capHir1', 'Broad/equCab2', 'CerSimSim1.0/cerSim1', 'ICGSC Felis_catus 6.2/felCat5',
            'Broad CanFam3.1/canFam3', 'MusPutFur1.0/musFur1', 'BGI-Shenzhen 1.0/ailMel1', 'Oros_1.0/odoRosDiv1',
            'LepWed1.0/lepWed1', 'ASM32557v1/pteAle1', 'Broad/pteVam1', 'ASM32734v1/myoDav1',
            'Broad Institute Myoluc2.0/myoLuc2', 'EptFus1.0/eptFus1', 'EriEur2.0/eriEur2', 'Broad/sorAra2',
            'ConCri1.0/conCri1',
            # Afrotheria subset
            'Broad/loxAfr3', 'EleEdw1.0/eleEdw1', 'Broad v1.0/triMan1', 'ChrAsi1.0/chrAsi1', 'Broad/echTel2',
            'OryAfe1.0/oryAfe1',
            # Mammal subset
            'Baylor/dasNov3', 'Broad/monDom5', 'WTSI Devil_ref v7.0/sarHar1', 'TWGS Meug_1.1/macEug2',
            'WUGSC 5.0.1/ornAna1',
            # Aves subset
            'F_cherrug_v1.0/falChe1', 'F_peregrinus_v1.0/falPer1', 'FicAlb1.5/ficAlb2', 'ASM38545v1/zonAlb1',
            'GeoFor_1.0/geoFor1', 'WashU taeGut324/taeGut2', 'PseHum1.0/pseHum1', 'WUSTL v6.3/melUnd1', 'AV1/amaVit1',
            'SMACv1.1/araMac1', 'Cliv_1.0/colLiv1', 'BGI_duck_1.0/anaPla1', 'ICGSC Gallus_gallus-4.0/galGal4',
            'TGC Turkey_2.01/melGal1',
            # Sarcopterygii subset
            'allMis0.2/allMis1', 'CheMyd_1.0/cheMyd1', 'v3.0.1/chrPic1', 'PelSin_1.0/pelSin1', 'ASM38561v1/apaSpi1',
            'Broad AnoCar2.0/anoCar2', 'JGI 7.0/xenTro7', 'Broad/latCha1',
            # Fish subset
            'Genoscope 8.0/tetNig2', 'FUGU5/fr3', 'version 1 of Takifugu flavidus genome/takFla1',
            'Broad oreNil1.1/oreNil2', 'NeoBri1.0/neoBri1', 'AstBur1.0/hapBur1', 'MetZeb1.1/mayZeb1',
            'PunNye1.0/punNye1',
            'NIG/UT MEDAKA1/oryLat2', 'Xiphophorus_maculatus-4.4.2/xipMac1', 'Broad/gasAcu1',
            'Genofisk GadMor_May2010/gadMor1', 'Zv9/danRer7', 'Astyanax_mexicanus-1.0.2/astMex1', 'LepOcu1/lepOcu1',
            'WUGSC 7.0/petMar2']
SPCS = [a.split('/')[-1] for a in SPCS]  # remove "full" species name if there was one

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

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N", "-": "-",
              "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}


####################################################################################################
# PROCESS PROTEIN INFORMATION
####################################################################################################

def trans_to_protein(fasta_file):
    """
    :param fasta_file: full path to a FASTA formatted file with chromosome, gene, protein, and transcript
                       information available in the header
    :return: dictionary of transcript ID -> (chromosome, gene, protein) IDs
    """

    trans_to_protloc = {}

    fasta_handle = gzip.open(fasta_file) if fasta_file.endswith('gz') else open(fasta_file)
    for fasta_line in fasta_handle:
        if fasta_line.startswith('>'):
            # remove the version for the protein, gene, and transcript Ensembl IDs:
            prot_id = fasta_line[fasta_line.find('prot:') + 5:-1].split()[0]  # .split('.')[0]
            gene_id = fasta_line[fasta_line.find('gene:') + 5:-1].split()[0]  # .split('.')[0]
            tran_id = fasta_line[fasta_line.find('transcript:') + 11:-1].split()[0].split('.')[0]
            chromosome = fasta_line[fasta_line.find('chromosome:' + GENOME_BUILD + ':'):].split()[0].split(':')[2]

            if tran_id not in trans_to_protloc:
                trans_to_protloc[tran_id] = set()
            trans_to_protloc[tran_id].add((chromosome, gene_id, prot_id))
    fasta_handle.close()

    return trans_to_protloc


########################################################################################################

def parse_alignment(fasta_file, alignment_file, output_dir):
    """
    :param alignment_file: full path to a FASTA formatted multiple alignment file BY PROTEIN, downloaded from:
                           http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/alignments/knownGene.exonAA.fa.gz
    :param output_dir: full path to a directory to store chromosome -> gene ID -> protein alignments and JSD values
    :return: None, but print success messages
    """

    # retrieve the ID mappings:
    tranid_to_protloc = trans_to_protein(fasta_file)

    # write out to files:
    current_files = None  # keep track of the file name we are writing to
    current_file_handles = None  # and the handle itself

    aln_handle = gzip.open(alignment_file) if alignment_file.endswith('gz') else open(alignment_file)
    for aln_line in aln_handle:
        if aln_line.startswith('>') and ('_' + BUILD_ALT_ID + '_') in aln_line:

            # close the previous file that we had been writing to
            if current_files:
                for handle in current_file_handles:
                    handle.close()

            # try to find the new file to write to:
            current_id = aln_line[1:-1].split('_')[0].split('.')[0]  # remove transcript version
            if current_id in tranid_to_protloc:
                current_files = []
                current_file_handles = []
                for chrom, gene_id, prot_id in tranid_to_protloc[current_id]:
                    if chrom not in os.listdir(output_dir):
                        call(['mkdir', output_dir + chrom])
                    if gene_id not in os.listdir(output_dir + chrom):
                        call(['mkdir', output_dir + chrom + '/' + gene_id])
                    current_files.append(output_dir + chrom + '/' + gene_id + '/' + prot_id + '.100way-alignment.fa.gz')

                    current_file_handles.append(
                        gzip.open(current_files[-1], 'a' if os.path.isfile(current_files[-1]) else 'w'))
            else:
                sys.stderr.write('No such ID: ' + current_id + '\n')
                current_files = None
                current_file_handles = None

        if current_file_handles:
            for handle in current_file_handles:
                handle.write(aln_line)
    aln_handle.close()

    # close last few remaining files if necessary:
    if current_files:
        for handle in current_file_handles:
            handle.close()

    # print success message
    sys.stderr.write('Wrote to ' + output_dir + '\n')


####################################################################################################
# PROCESS EXON INFORMATION
####################################################################################################

def get_exon_locs(protfile,
                  chromosome_set=tuple(['chr' + chr_val for chr_val in map(str, range(1, 23)) + ['X', 'Y']])):
    """
    :param protfile: full path to a fasta-formatted file containing protein sequences where the headers
                     of the sequences have been updated to include exon ranges in order
    :param chromosome_set: set of chromosome names to consider
    :return: an IntervalTree of exon locations on the appropriate chromosome
    """

    exon_locations = {chrom_id: IntervalTree() for chrom_id in chromosome_set}

    if not os.path.isfile(protfile):
        sys.stderr.write('Could not open ' + protfile + '\n')
        return exon_locations

    exon_intervals_by_gene = {}

    fasta_handle = gzip.open(protfile) if protfile.endswith('gz') else open(protfile)
    for seq_line in fasta_handle:
        if seq_line.startswith('>'):
            loc = seq_line[seq_line.find(GENOME_BUILD+':')+len(GENOME_BUILD+':'):].split()[0]
            chrom_id = loc.split(':')[0]
            chrom_id = ('' if chrom_id.startswith('chr') else 'chr') + chrom_id

            if chrom_id not in chromosome_set:
                continue

            exons = loc.split(':')[1].replace('join(', '').replace('complement(', '').replace(')', '').split(',')
            exons = [tuple(sorted([int(exon_range.split('..')[0]), int(exon_range.split('..')[1])]))
                     for exon_range in exons]

            gene_id = seq_line[seq_line.find('gene:') + 5:].split()[0]  # e.g., ENSG00000000457.9

            if chrom_id not in exon_intervals_by_gene:
                exon_intervals_by_gene[chrom_id] = {}
            if gene_id not in exon_intervals_by_gene[chrom_id]:
                exon_intervals_by_gene[chrom_id][gene_id] = set()

            # Add each POSITIVE exon SEPARATELY!
            for exon_start, exon_end in exons:  # these are 1-indexed, but we need 0-indexed...
                if exon_start > -1 and exon_end > -1:
                    exon_intervals_by_gene[chrom_id][gene_id].add((exon_start, exon_end))
    fasta_handle.close()

    # Store all the information as interval trees:
    for chrom_id in exon_intervals_by_gene.keys():
        for gene_id in exon_intervals_by_gene[chrom_id].keys():
            for exon_start, exon_end in exon_intervals_by_gene[chrom_id][gene_id]:
                exon_locations[chrom_id].addi(exon_start, exon_end + (1 if exon_end == exon_start else 0),
                                              gene_id + '_' + str(exon_start) + '-' + str(exon_end))

    return exon_locations


####################################################################################################

def concatenate_maf_to_fasta(input_maf, output_fasta, chrom_start_pos, chrom_end_pos, species=SPCS):
    """
    :param input_maf: full path to a MultiZWay alignment file (that we just created for a single block)
    :param output_fasta: full path to resulting fasta multiple alignment file
    :param chrom_start_pos: starting index in the chromosome for this particular species for this alignment block
    :param chrom_end_pos: ending index in the chromosome for this particular species for this alignment block
    :param species: all the species that should be represented in file
    :return: None
    """

    fill = ""  # how to separate blocks?
    wrap = -1  # maximum line length in alignment (or no maximum)

    fasta_handle = gzip.open(output_fasta, 'w') if output_fasta.endswith('gz') else open(output_fasta, 'w')

    texts = {current_species: [] for current_species in species}  # store the complete sequence for each species

    # Read in the maf file !
    maf_handle = gzip.open(input_maf) if input_maf.endswith('gz') else open(input_maf)
    for maf_aln in AlignIO.parse(maf_handle, "maf"):  # for each alignment block
        observed_species, sequence_lengths = set(), set()
        for seq_rec in maf_aln:
            current_species = seq_rec.id.split('.')[0]
            if current_species in species:
                texts[current_species].append(str(seq_rec.seq))
                observed_species.add(current_species)
                sequence_lengths.add(len(str(seq_rec.seq)))
        for current_species in species:
            if current_species not in observed_species:
                texts[current_species].append("-" * max(sequence_lengths))
    maf_handle.close()

    if False:
        maf_reader = maf.Reader(maf_handle)
        for maf_line in maf_reader:
            for current_species in species:
                sequence_component = maf_line.get_component_by_src_start(current_species)
                if sequence_component:
                    texts[current_species].append(sequence_component.text)
                else:
                    texts[current_species].append("-" * maf_line.text_size)
        maf_handle.close()

    # Finally, print out the concatenated sequence by species in order
    for current_species in species:
        fasta_handle.write('>' + current_species + (' ' + str(chrom_start_pos + 1) + ':' + str(chrom_end_pos + 1)
                                                    if current_species == BUILD_ALT_ID else '') + '\n')
        complete_matching_sequence = fill.join(texts[current_species])
        if wrap <= 0:
            fasta_handle.write(complete_matching_sequence + '\n')

        else:
            p = 0  # starting index of entire alignment for each line
            while p < len(complete_matching_sequence):
                fasta_handle.write(complete_matching_sequence[p:min(p + wrap, len(complete_matching_sequence))] + '\n')
                p += wrap

    fasta_handle.close()
    sys.stderr.write('Reformatted '+input_maf+' to '+output_fasta+'\n')


####################################################################################################

def extract_exon_fasta_from_maf(chroms, input_directory, output_directory):
    """
    :param chroms: set of chromosomes to run on
    :return: None, but print success message upon successful rewrite of fasta-formatted cDNA alignments per EXON
    """

    ensembl_to_ucscgb_chr_mapping = {'chrKI270711.1': 'chr1_KI270711v1_random',
                                     'chrKI270713.1': 'chr1_KI270713v1_random',
                                     'chrKI270721.1': 'chr11_KI270721v1_random',
                                     'chrGL000009.2': 'chr14_GL000009v2_random',
                                     'chrGL000194.1': 'chr14_GL000194v1_random',
                                     'chrKI270726.1': 'chr14_KI270726v1_random',
                                     'chrKI270727.1': 'chr15_KI270727v1_random',
                                     'chrKI270728.1': 'chr16_KI270728v1_random',
                                     'chrGL000205.2': 'chr17_GL000205v2_random',
                                     'chrKI270731.1': 'chr22_KI270731v1_random',
                                     'chrKI270734.1': 'chr22_KI270734v1_random',
                                     'chrGL000195.1': 'chrUn_GL000195v1',
                                     'chrGL000213.1': 'chrUn_GL000213v1',
                                     'chrGL000218.1': 'chrUn_GL000218v1',
                                     'chrGL000219.1': 'chrUn_GL000219v1',
                                     'chrMT': 'chrM'}

    for chrom_id in chroms:
        sys.stderr.write('Getting exon locations for genes on chromosome(s) ' + ','.join([chrom_id]) + '\n')
        # All gene locations on this chromosome
        gene_locs = get_exon_locs(data_path+'ensembl/Homo_sapiens.'+GENOME_BUILD+'/Homo_sapiens.'+GENOME_BUILD +
                                  '.pep.all.withgenelocs.verified.fa.gz', [chrom_id])
        sys.stderr.write('Finished!\n')

        full_chrom_name = ensembl_to_ucscgb_chr_mapping.get(chrom_id, chrom_id)

        maf_file = input_directory + full_chrom_name + '.maf.gz'
        gene_path = output_directory + chrom_id.replace('chr', '') + '/'

        if not os.path.isfile(maf_file):
            sys.stderr.write('Could not open ' + maf_file + '\n')
            continue

        if not os.path.isdir(gene_path):
            sys.stderr.write('No such directory ' + gene_path + '\n')
            continue
            # os.system('mkdir ' + gene_path)

        newfiles = {}  # The new smaller output .maf files and the regions they cover

        current_block = []  # Temporary list to keep track of each alignment chunk

        exon_handle = gzip.open(maf_file) if maf_file.endswith('gz') else open(maf_file)

        maf_header = []  # Each .maf file must start with ##maf ...
        for exon_line in exon_handle:
            if not exon_line.startswith('a score'):
                maf_header.append(exon_line)
            else:
                current_block.append(exon_line)
                break

        sys.stderr.write('Processed ' + maf_file + ' header...\n')

        for exon_line in exon_handle:

            # Write out the current block:
            if exon_line.startswith('\n') and len(current_block) > 0:
                current_block.append(exon_line)
                startloc, seqlen = map(int, [a for a in current_block
                                             if BUILD_ALT_ID + '.' + full_chrom_name in a][0].split()[2:4])

                # Remember that the .maf files are 0-indexed whereas the exon locations are 1-indexed.
                for start, end, exon_id in gene_locs[chrom_id][startloc + 1:startloc + seqlen + 1]:

                    gene_id = exon_id.split('_')[0]

                    genefile = gene_path + gene_id + '/' + exon_id + '.alignment.maf'

                    if genefile not in newfiles:
                        out_handle = open(genefile, 'w')
                        map(out_handle.write, maf_header)
                        newfiles[genefile] = [startloc, startloc + seqlen]
                    else:
                        out_handle = open(genefile, 'a')
                        newfiles[genefile][1] = startloc + seqlen

                    map(out_handle.write, current_block)
                    out_handle.close()
                current_block = []

            else:
                current_block.append(exon_line)
        exon_handle.close()

        # Finally, convert to a fasta file !
        for exon_maf_file in sorted(list(newfiles.keys())):
            sys.stderr.write('Concatenating ' + exon_maf_file + ' to fasta...\n')
            concatenate_maf_to_fasta(exon_maf_file,
                                     exon_maf_file.replace('.alignment.maf', '.aln.fa.gz'),
                                     newfiles[exon_maf_file][0],
                                     newfiles[exon_maf_file][1])
            os.system('rm ' + exon_maf_file)


####################################################################################################

def run_per_exon(chroms, function, aln_directory):
    """
    :param chroms: set of chromosomes to run on
    :param function: pointer to a function to call in the innermost loop (per protein)
    :param aln_directory: full path to a directory where alignment files have been stored by chromosome
    :return: None, but print appropriate success message
    """

    total_failed = 0
    total_processed = 0
    for chrom_id in chroms:
        if not os.path.isdir(aln_directory + chrom_id.replace('chr', '')):
            continue

        gene_list = os.listdir(aln_directory + chrom_id.replace('chr', '') + '/')

        for gene_id in gene_list:
            if not os.path.isdir(aln_directory + chrom_id.replace('chr', '') + '/' + gene_id + '/'):
                continue

            for prot_id in list(set([a.replace('.exons.txt', '').replace('.100way-alignment.fa.gz', '') for a in
                            os.listdir(aln_directory + chrom_id.replace('chr', '') + '/' + gene_id)
                            if a.endswith('.exons.txt') or a.endswith('.100way-alignment.fa.gz')])):

                # full path to the exons file (sorted exons)
                exonfile = aln_directory + chrom_id.replace('chr', '') + '/' + gene_id + '/' + prot_id + '.exons.txt'
                if not os.path.isfile(exonfile) and function in ['stitch_exons', 'repair_cdna_prots']:
                    continue

                if function == 'calc_jsd_conservation':
                    if os.path.isfile(exonfile):
                        with open(exonfile) as exon_handle:
                            header = exon_handle.next()
                    else:
                        header = '>'+prot_id+'\n'

                # set of all exon alignments (we'll see which ones we need for this particular protein)
                if function == 'stitch_exons':
                    alignments = {(a.split('_')[1].split('-')[0],
                                   a.split('-')[1].split('.')[0]):
                                      aln_directory + chrom_id.replace('chr', '') + '/' + gene_id + '/' + a
                                  for a in os.listdir(aln_directory + chrom_id.replace('chr', '') + '/' + gene_id + '/')
                                  if a.startswith(gene_id) and a.endswith('.aln.fa.gz')}

                # input protein and cdna files
                protein_file = aln_directory + chrom_id.replace('chr', '') + '/' + gene_id + '/' + prot_id + '.prot.fa'
                cdna_file = aln_directory + chrom_id.replace('chr', '') + '/' + gene_id + '/' + prot_id + '.cdna.fa'

                # output alignment file
                alignment_file = aln_directory + chrom_id.replace('chr', '') + '/' + gene_id + '/' + \
                                 prot_id + '.100way-alignment.fa.gz'

                # output JSD file:
                jsd_file = aln_directory + chrom_id.replace('chr', '') + '/' + gene_id + '/' + prot_id + '.jsd.txt'

                if function == 'stitch_exons':
                    stitch_exons(exonfile, alignments, alignment_file)
                elif function == 'calc_jsd_conservation':
                    calc_jsd_conservation(alignment_file, jsd_file, header)
                elif function == 'repair_cdna_prots':
                    fix_cdna_prots(exonfile, cdna_file, protein_file)
                elif function == 'check_alignments':
                    passing = check_alignments(alignment_file, protein_file)
                    total_processed += 1
                    if not passing:
                        total_failed += 1

    print str(total_failed) + '/' + str(total_processed)


####################################################################################################

def create_protein_alignment(chroms, aln_directory):
    """
    :param chroms: set of chromosomes to run on
    :param aln_directory: full path to a directory where alignment files have been stored by chromosome
    :return: None, but print success message upon successful write of fasta-formatted cDNA alignments per protein
    """

    run_per_exon(chroms, 'stitch_exons', aln_directory)


####################################################################################################

def translate_rna(rna_seq):
    """
    :param rna_seq: input cDNA sequence string
    :return: translation of the cDNA sequence to protein using GENCODE; unknown codons or accidental
             stop codons are translated as 'X'
    """

    prot_seq = []
    for i in range(0, len(rna_seq), 3):
        codon = rna_seq[i:i + 3]
        if len(codon) == 3:
            if codon not in GENCODE and '-' in codon:
                prot_seq.append('-')
            elif codon not in GENCODE and 'N' in codon:
                prot_seq.append('X')
            elif GENCODE[codon] == "_" and i < len(rna_seq) - 5:  # stop codon somewhere in the middle..
                prot_seq.append('X')
            else:
                prot_seq.append(GENCODE[codon])
    return ''.join(prot_seq)


####################################################################################################

def stitch_exons(exonfile, exon_alignment_files, aln_outfile, error_check=False):
    """
    :param exonfile: full path to the exons file (sorted exons)
    :param exon_alignment_files: set of all exon alignments (we'll see which ones we need for this particular protein)
    :param aln_outfile: full path to where to write a multiple protein alignment
    :param error_check: whether to print intermediate steps (for debugging purposes)
    :return: None, but given the exon files we have already created and the fasta formatted maf blocks from
             whole genome alignments, create an amino acid multiple sequence alignment for the protein of interest
    """

    sys.stderr.write('Attempting to run on ' + exonfile + '...\n')

    cdna_per_spc = {spc: [] for spc in SPCS}  # species -> cDNA sequence (to be translated at the end)
    negative_strand = False  # are we looking at the complementary strand for the ENTIRE PROTEIN?

    exon_handle = gzip.open(exonfile) if exonfile.endswith('gz') else open(exonfile)
    for exon_line in exon_handle:
        if exon_line.startswith('>'):  # store whether this is a negative strand sequence
            if 'complement' in exon_line:
                negative_strand = True
            continue

        # for each exon (one per line):
        complement = False  # is this particular exon already reversed?
        exon_loc, exon_seq = exon_line[:-1].split('\t')

        exon_start, exon_end = map(int, exon_loc.split(':'))

        # swap the start and end if necessary
        if exon_end < exon_start or (exon_end < 0 and exon_end == exon_start and negative_strand):
            complement = True
            tmp = exon_start
            exon_start = exon_end
            exon_end = tmp

        # if we can't find the corresponding exon alignment, fluff this exon with gaps, and move on
        if (str(exon_start), str(exon_end)) not in exon_alignment_files:
            sys.stderr.write('No alignment file for ' + exonfile + ': ' + str(exon_start) + '-' + str(exon_end) + '\n')
            cdna_per_spc[BUILD_ALT_ID].append(''.join([COMPLEMENT[seq_index] for seq_index in list(exon_seq)])
                                              if complement else exon_seq)
            for spc in SPCS[1:]:
                cdna_per_spc[spc].append('-' * len(exon_seq))
            continue

        # otherwise, keep track of the exon information
        with gzip.open(exon_alignment_files[(str(exon_start), str(exon_end))]) as exon_aln_handle:
            human_block_location = exon_aln_handle.next()[:-1].split()[1]
            block_start, block_end = map(int, human_block_location.split(':'))  # this is ALWAYS in the correct order

            human_sequence = exon_aln_handle.next()[:-1]

            # these indices will correspond to the indices into the sequence with gaps (i.e. not continuous)
            nongap_indices = [(nt, seq_index) for (nt, seq_index) in
                              [(nt, seq_index) for seq_index, nt in enumerate(list(human_sequence))] if nt != '-']

            if error_check:
                sys.stderr.write('------------------\n')
                sys.stderr.write('exon_start = ' + str(exon_start) + '\n')
                sys.stderr.write('block_start = ' + str(block_start) + '\n')
                sys.stderr.write('exon_end = ' + str(exon_end) + '\n')
                sys.stderr.write('block_end = ' + str(block_end) + '\n')

            # if the block is SMALLER than our desired exon, pad out the block accordingly:
            fluff_left = 0 if exon_start >= block_start else block_start - exon_start + 1
            fluff_right = 0 if exon_end <= block_end else exon_end - block_end + 1

            if complement:
                tmp = fluff_left
                fluff_left = fluff_right
                fluff_right = tmp

            # shave off excess nucleotides from front and back:
            start_exon_index = 0 if exon_start <= block_start else exon_start - block_start  # shave off excess in front
            end_exon_index = len(nongap_indices) if exon_end >= block_end else len(nongap_indices) - (
                block_end - exon_end - 1)

            if error_check:
                sys.stderr.write('fluff_left: ' + str(fluff_left) + '\n')
                sys.stderr.write('fluff_right: ' + str(fluff_right) + '\n')
                sys.stderr.write('start_block_index: ' + str(start_exon_index) + '\n')
                sys.stderr.write('end_block_index: ' + str(end_exon_index) + '\n')

            # get ALL indices (not just ranges)
            this_exon = nongap_indices[start_exon_index:end_exon_index]
            indices = [index for nuc, index in this_exon]
            if complement:
                indices = indices[::-1]

            fixed = True  # have we "fixed" a problematic exon?
            if ''.join([human_sequence[seq_index] for seq_index in indices]).upper() != \
                exon_seq.upper()[fluff_left:len(exon_seq) - fluff_right]:

                if error_check:
                    sys.stderr.write(
                        'us-----: ' + ''.join([human_sequence[seq_index] for seq_index in indices]).upper() + '\n')
                    sys.stderr.write('correct: ' + exon_seq.upper()[fluff_left:len(exon_seq) - fluff_right] + '\n')

                fixed = False

                newexons = {
                    'remove_left': nongap_indices[start_exon_index + 1:end_exon_index],
                    'remove_right': nongap_indices[start_exon_index:end_exon_index - 1],
                    'add_left': nongap_indices[start_exon_index - (1 if start_exon_index > 0 else 0):end_exon_index],
                    'add_right': nongap_indices[start_exon_index:end_exon_index +
                                                                 (1 if end_exon_index < len(nongap_indices) else 0)],
                    'shift_left': nongap_indices[start_exon_index-(1 if start_exon_index > 0 else 0):end_exon_index-1],
                    'shift_right': nongap_indices[start_exon_index+1:end_exon_index +
                                                                     (1 if end_exon_index < len(nongap_indices) else 0)]
                }

                # for each TYPE of fix (see keys above) and the resulting new exon,
                for fix_type, new_exon in newexons.items():
                    fixed_indices = [index for nuc, index in new_exon]
                    if complement:
                        fixed_indices = fixed_indices[::-1]

                    # keep track of original "fluff" left/right values incase they are reset in the next step:
                    orig_fluff_left = fluff_left
                    orig_fluff_right = fluff_right

                    # try to add left?
                    if fix_type in {'add_left', 'shift_left'} and start_exon_index <= 0:
                        if complement:
                            fluff_right += 1
                        else:
                            fluff_left += 1

                    # try to add right?
                    if fix_type in {'add_right', 'shift_right'} and end_exon_index >= len(nongap_indices):
                        if complement:
                            fluff_left += 1
                        else:
                            fluff_right += 1

                    if error_check:
                        sys.stderr.write(fix_type.split()[0][0] + fix_type.split()[1][0] + '-----: ' +
                                         ''.join(
                                             [human_sequence[seq_index] for seq_index in fixed_indices]).upper() + '\n')

                    if ''.join([human_sequence[seq_index] for seq_index in fixed_indices]).upper() == \
                        exon_seq.upper()[fluff_left:len(exon_seq) - fluff_right]:
                        indices = fixed_indices
                        fixed = True
                        final_fix_type = fix_type if not complement else (
                            fix_type.replace('left', 'right') if 'left' in fix_type else
                            fix_type.replace('right', 'left'))
                        sys.stderr.write('Fixed by ' + final_fix_type + '\n')
                        break

                    fluff_left = orig_fluff_left
                    fluff_right = orig_fluff_right

            if not fixed:
                sys.stderr.write('Exon ' + exonfile + ': ' + str(exon_start) + '-' + str(exon_end) +
                                 ' did not match:\n' +
                                 ''.join([human_sequence[seq_index] for seq_index in indices]).upper() + '\n' +
                                 exon_seq.upper() + '\n')

                for fix_type, new_exon in newexons.items():
                    final_fix_type = fix_type if not complement else (
                        fix_type.replace('left', 'right') if 'left' in fix_type else
                        fix_type.replace('right', 'left'))
                    sys.stderr.write(final_fix_type + '\n' +
                                     ''.join([human_sequence[seq_index] for seq_index in
                                              [exon_index[1] for exon_index in new_exon][::-1
                                              if complement else 1]]).upper() +
                                     '\n#\n')

                cdna_per_spc[BUILD_ALT_ID].append(
                    ''.join([COMPLEMENT[seq_index] for seq_index in list(exon_seq)]) if complement else exon_seq)
                for spc in SPCS[1:]:
                    cdna_per_spc[spc].append('-' * len(exon_seq))
                continue

            cdna_per_spc[BUILD_ALT_ID].append(('-' * fluff_left) +
                                              ''.join([COMPLEMENT[human_sequence[seq_index]]
                                                       for seq_index in indices]) if complement else
                                              ''.join([human_sequence[seq_index] for seq_index in indices]) +
                                              ('-' * fluff_right))

            # for the remaining lines in the exon alignment (we already processed the first two lines = human)
            for aln_line in exon_aln_handle:
                if aln_line.startswith('>'):
                    currspc = aln_line[1:-1]
                    currseq = exon_aln_handle.next()[:-1]
                    cdna_per_spc[currspc].append(('-' * fluff_left) +
                                                 ''.join(COMPLEMENT[currseq[i]] for i in indices) if complement else
                                                 ''.join(currseq[i] for i in indices) +
                                                 ('-' * fluff_right))
    exon_handle.close()

    out_handle = gzip.open(aln_outfile, 'w') if aln_outfile.endswith('gz') else open(aln_outfile, 'w')
    for spc in SPCS:
        out_handle.write('>' + spc + '\n' + translate_rna(''.join(cdna_per_spc[spc]).upper()) + '\n')
    out_handle.close()
    sys.stderr.write('Wrote multiple protein alignment to ' + aln_outfile + '\n')


####################################################################################################

def fix_cdna_prots(exonfile, cdnafile, protfile):
    """
    :param exonfile: full path to an exons .txt file
    :param cdnafile: full path to a FASTA-formatted cDNA file
    :param protfile: full path to a FASTA-formatted protein file
    :return: none, but "fix" the cDNA and protein files in some cases where the complement was not
             taken properly. You will likely never need to call this function
    """

    rnaseq = []
    with open(exonfile) as exon_handle:
        header = exon_handle.next()
        complement = 'complement' in header

        for exon_line in exon_handle:
            exon = exon_line[:-1].split()[1]
            rnaseq.append(''.join([COMPLEMENT[i] for i in list(exon)]) if complement else exon)

    cdna_handle = open(cdnafile, 'w')
    cdna_handle.write(header + ''.join(rnaseq) + '\n')
    cdna_handle.close()

    prot_handle = open(protfile, 'w')
    prot_handle.write(header + translate_rna(''.join(rnaseq).upper())[:-1] + '\n')
    prot_handle.close()


####################################################################################################
# CALCULATE CONSERVATION FROM PROTEIN ALIGNMENTS
####################################################################################################

def calculate_jsd_scores(chroms, aln_directory):
    """
    :param chroms: set of chromosomes to run on
    :param aln_directory: full path to a directory where alignment files have been stored by chromosome
    :return: None, but print success message upon successful write of fasta-formatted cDNA alignments per protein
    """

    run_per_exon(chroms, 'calc_jsd_conservation', aln_directory)


####################################################################################################

def KLdiv(po, qo):
    """
    :param po: vector of values
    :param qo: vector of values (should be same length as po)
    :return: the Kullback-Leibler divergence of two vectors, KL(p||q)
    """

    assert len(po) == len(qo), \
        "Cannot compute KL divergence of vectors of unequal length."
    p = [float(i) for i in po]
    q = [float(i) for i in qo]
    return sum([p[x] * math.log(p[x] / q[x], 2) for x in range(len(p)) \
                if p[x] != 0. and p[x] != 0 and q[x] != 0. and q[x] != 0])


####################################################################################################

def JSdiv(po, qo, wt=0.5):
    """
    :param po: vector of values
    :param qo: vector of values (should be same length as po)
    :param wt: weight of the first vector
    :return: the Jensen-Shannon divergence of two vectors, JS(p||q)
    """

    assert len(po) == len(qo), \
        "Cannot compute JS divergence of vectors of unequal length."
    p = [float(i) for i in po]
    q = [float(i) for i in qo]

    # Take weighted average of KL divergence
    av = [wt * p[x] + (1 - wt) * q[x] for x in xrange(len(p))]
    return wt * KLdiv(p, av) + (1 - wt) * KLdiv(q, av)


####################################################################################################

def JSdiv_column(residues, rand=False):
    """Compute the Jensen-Shannon divergence of vector r with background
    frequency. Default is Blosum62 background frequencies, but can use
    random if specified"""

    aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
          'T', 'W', 'Y', 'V']  # 20 amino acids, not including gaps...
    residues = [r for r in residues if r in aa]  # Remove gaps (-) and unknowns (X)

    if len(residues) < 2: return 0  # Too few seqs in this particular column

    # Background frequencies of amino acids (random or BLOSUM62 matrix):
    q = [1. / len(aa)] * len(aa) if rand else \
        [0.074, 0.052, 0.045, 0.054, 0.025, 0.034, 0.054, 0.074,
         0.026, 0.068, 0.099, 0.058, 0.025, 0.047, 0.039, 0.057,
         0.051, 0.013, 0.032, 0.073]

    fqs = [0.00001 + float(residues.count(s)) / len(residues) for s in aa]
    freqs = [f / sum(fqs) for f in fqs]

    assert str(sum(q)) == '1.0' and str(sum(freqs)) == '1.0', \
        "Prob. vectors do not sum to 1"

    return JSdiv(freqs, q)


####################################################################################################

def calc_jsd_conservation(alnfile, outfile, header=''):
    """
    :param alnfile: full path to a protein multiple sequence alignment
    :param outfile: full path to a file to write the Jensen-Shannon divergences from a Blosum 62 background (as a
                    measure of conservation) of each column in the alignment
    :param header: protein Ensembl description header (obtained from the exon file...
    :return: None, but write success message
    """

    # 100-vertebrate protein alignment with human protein sequence:
    if not os.path.isfile(alnfile):
        return

    # Store the full alignment across all species (to calculate column-based scores)
    seqs = {}
    mainspc = ''
    aln_handle = gzip.open(alnfile) if alnfile.endswith('gz') else open(alnfile)
    for aln_line in aln_handle:
        if aln_line.startswith('>'):
            try:
                spc = aln_line[1:-1].split()[0].split('_')[1]  # species name (e.g., hg38...)
            except IndexError:
                spc = aln_line[1:-1].split()[0]  # if the alignments are not taken directly from UCSC GB
            if spc not in seqs:
                seqs[spc] = []
            if mainspc == '':
                mainspc = spc
            seqs[spc].append(aln_handle.next().strip())  # full protein sequence (NOTE: ends with stop codon)
    aln_handle.close()

    for spc in seqs.keys():
        seqs[spc] = ''.join(seqs[spc])
        if seqs[spc].endswith('Z'):
            seqs[spc] = seqs[spc][:-1]  # remove "stop" codon

    # NOW, calculate the column-based JSD...
    aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
          'T', 'W', 'Y', 'V']  # 20 amino acids, not including gaps...

    # Output file containing the JSD conservation values for each position
    out_handle = open(outfile, 'w')
    out_handle.write(header)

    for col_index in xrange(len(seqs[mainspc])):

        # residue value at this position in the human sequence
        curr_aa = seqs[mainspc][col_index]

        # full ungapped column
        curr_col = [curr_seq[col_index] for curr_seq in seqs.values()
                    if col_index < len(curr_seq) and curr_seq[col_index] in aa]

        # overall conservation score
        curr_val = len(curr_col) / 100. * JSdiv_column(curr_col)
        out_handle.write('\t'.join([str(col_index), str(curr_aa), str(curr_val)]) + '\n')

    out_handle.close()
    sys.stderr.write('Wrote ' + outfile + '\n')


####################################################################################################

def check_alignments(alignment_file, protein_file):
    """
    :param chrom: chromosome ID
    :param gene: gene ID
    :param prot: protein ID
    :return: whether the alignment sequence matches (well enough) the expected protein sequence
    """

    with gzip.open(alignment_file) as aln_file:
        aln_file.next()
        aln_seq = aln_file.next().strip()

    with open(protein_file) as prot_file:
        prot_file.next()
        prot_seq = prot_file.next().strip()

    distance_allowance = 0.05 * len(prot_seq)
    hamming_distance = sum([1 if prot_seq[i] != (aln_seq[i] if i < len(aln_seq) else '-') else 0
                            for i in xrange(len(prot_seq))])

    if hamming_distance > distance_allowance:
        sys.stderr.write('hamming distance > '+str(distance_allowance)+': '+protein_file+'\n')
        return False
    else:
        return True


####################################################################################################

def create_100way_jsd_files(domainfile, weightfile, exon_dir):
    """
    :param domainfile: full path to a file to write OUT to (track locations per protein)
    :param weightfile: full path to a file to write OUT ot (weights per track position by track name)
    :param exon_dir: full path to a directory containing subdirectories by chromosome and JSD values
    :return: None, but print success message upon successful write of two output files
    """

    weightfile_handle = gzip.open(weightfile, 'w') if weightfile.endswith('gz') else open(weightfile, 'w')
    weightfile_handle.write('\n'.join(
        ['# Continuous conservation weights across entire protein sequences, calculated using UCSC Genome Browser\'s',
         '#   100-vertebrate multiple alignment (http://hgdownload.soe.ucsc.edu/goldenPath/' + BUILD_ALT_ID +
         '/multiz100way/), downloaded December 3, 2018.',
         '# Nucleotide alignments for all protein isoforms were extracted and translated to create per-protein ' +
         'multiple alignments.',
         '# Conservation scores per column were calculated as the Jensen-Shannon divergence of the non-gapped residues',
         '#   against a Blosum 62 background (https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt), multiplied',
         '#   by the fraction of non-gapped residues in the column.',
         '\t'.join(['#EnsemblProtID_JSD', 'LigandTypeAlwaysCONS_', '1-Indexed-ProteinPosition', 'ConservationScore',
                    'NumInstancesAlways10', '1-Indexed-ProteinPositionAndResidue']) + '\n']))

    domainfile_handle = gzip.open(domainfile, 'w') if domainfile.endswith('gz') else open(domainfile, 'w')
    domainfile_handle.write(
        '\n'.join(['# All proteins for which conservation scores could be determined, using the UCSC Genome Browser\'s',
                   '#   100-vertebrate multiple alignment (http://hgdownload.soe.ucsc.edu/goldenPath/' + BUILD_ALT_ID +
                   '/multiz100way/), downloaded December 3, 2018.',
                   '\t'.join(['#EnsemblProtID', 'DomainID', 'matchstate:AA-0-index:AA-value\n'])]))

    for chrom in sorted(os.listdir(exon_dir)):
        if not os.path.isdir(exon_dir + chrom):
            continue

        for gene_id in sorted(os.listdir(exon_dir + chrom)):
            if not os.path.isdir(exon_dir + chrom + '/' + gene_id):
                continue

            for prot_id in sorted([a.replace('.jsd.txt', '') for a in os.listdir(exon_dir + chrom + '/' + gene_id)
                                   if a.endswith('.jsd.txt')]):
                conservation_score_available = False
                fullprotein = []

                exon_handle = open(exon_dir + chrom + '/' + gene_id + '/' + prot_id + '.jsd.txt')
                for exon_line in exon_handle:
                    if exon_line.startswith('>'):
                        continue
                    protpos, residue, cons = exon_line.strip().split('\t')[:3]
                    if len(protpos) < 1 or len(residue) < 1 or len(cons) < 1:
                        continue

                    cons = max(0., float(cons))
                    fullprotein.append(str(int(protpos) + 1) + ':' + protpos + ':' + residue)

                    if cons > 0:
                        conservation_score_available = True
                        weightfile_handle.write(
                            '\t'.join([prot_id + '_100way_JSD', '100WAY_CONS_', str(int(protpos) + 1), str(cons),
                                       '10', str(int(protpos) + 1) + ':' + residue]) + '\n')
                exon_handle.close()

                if conservation_score_available:
                    # last "amino acid" is SOMETIMES a stop codon, so remove it to prevent failures
                    if fullprotein[-1].split(':')[2] == '_': fullprotein = fullprotein[:-1]
                    domainfile_handle.write('\t'.join([prot_id, prot_id + '_100way_JSD', ','.join(fullprotein)]) + '\n')

    domainfile_handle.close()
    weightfile_handle.close()


####################################################################################################
# MAIN
####################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Convert multiple alignment into separate gene-based files.')
    parser.add_argument('--maf_directory', type=str, default=data_path+'ucscgb/'+BUILD_ALT_ID+'alignment/',
                        help='Full path to a directory where .maf files should be stored / can be found')
    parser.add_argument('--aln_directory', type=str, default=data_path+'ensembl/Homo_sapiens.'+GENOME_BUILD+'/exons/',
                        help='Full path to a directory where 100way alignment files have been saved')
    parser.add_argument('--outdir', type=str, default=data_path+'ucscgb/'+BUILD_ALT_ID+'alignment/',
                        help='Full path to where two output track files should be written.')
    parser.add_argument('--chromosome', type=str,
                        help='Chromosome number in ' + GENOME_BUILD + ' (i.e., 1-22, X, Y).',
                        default='chr1')

    parser.add_argument('--download_mafs', dest='download_mafs', action='store_true', default=False,
                        help='Download MultiZ-100way alignments, per-chromosome, from the UCSC Genome Browser')
    parser.add_argument('--download_fasta', dest='download_fasta', action='store_true', default=False,
                        help='Download all known gene protein alignments from the UCSC Genome Browser')
    parser.add_argument('--extract_protein_alignments', dest='extract_protein_alignments', action='store_true',
                        default=False, help='Extract protein alignments directly from fasta UCSC alignment file.')
    parser.add_argument('--create_exon_alignments', dest='create_exon_alignments', action='store_true', default=False,
                        help='Create exon alignments in FASTA format from the .maf blocks')
    parser.add_argument('--create_protein_alignments', dest='create_protein_alignments', action='store_true',
                        default=False, help='Create protein alignments from the DNA exon alignments in the last step.')
    parser.add_argument('--compute_jsd', dest='compute_jsd', action='store_true', default=False,
                        help='Use protein alignments to compute the per-column JSD conservation values.')
    parser.add_argument('--create_track_files', dest='create_track_files', action='store_true', default=False,
                        help='Create two input files required for PertInInt describing the tracks and ' +
                             'per-position weights corresponding to conservation in those tracks.')
    args = parser.parse_args()

    # create directories as needed:
    for new_dir in ['ensembl', 'ensembl/Homo_sapiens.' + GENOME_BUILD, 'ucscgb',
                    'ucscgb/' + BUILD_ALT_ID + 'alignment/', 'ucscgb/' + BUILD_ALT_ID + 'alignment/mafs/']:
        if not os.path.isdir(data_path + new_dir):
            call(['mkdir', data_path + new_dir])

    # ----------------------------------------------------------------------------------------------------
    if args.download_fasta:
        """Download the fasta-formatted multiple protein alignments from the UCSC Genome Browser"""
        if not os.path.isdir(args.maf_directory):
            call(['mkdir', args.maf_directory])
        if not args.maf_directory.endswith('/'):
            args.maf_directory += '/'

        sys.stderr.write('Downloading knownGene.exonAA.fa.gz... ')
        call(['wget',
              'http://hgdownload.soe.ucsc.edu/goldenPath/'+BUILD_ALT_ID+'/multiz100way/alignments/knownGene.exonAA.fa.gz',
              '-O', args.maf_directory+'knownGene.exonAA.fa.gz'])

    # ----------------------------------------------------------------------------------------------------
    elif args.extract_protein_alignments:
        """Extract protein alignments from the protein alignment file downloaded from the UCSC genome browser"""
        parse_alignment(data_path+'ensembl/Homo_sapiens.'+GENOME_BUILD+'/Homo_sapiens.'+GENOME_BUILD +
                                  '.pep.all.withgenelocs.verified.fa.gz',
                        args.maf_directory+'knownGene.exonAA.fa.gz',
                        args.aln_directory)

    # ----------------------------------------------------------------------------------------------------
    elif args.download_mafs:
        """Download all MultiZ alignment files from the UCSC Genome Browser"""
        if not os.path.isdir(args.maf_directory):
            call(['mkdir', args.maf_directory])
        if not args.maf_directory.endswith('/'):
            args.maf_directory += '/'
        if not os.path.isdir(args.maf_directory+'mafs/'):
            call(['mkdir', args.maf_directory+'mafs/'])
        if not os.path.isfile(args.maf_directory+'directory_list.html'):
            call(['wget', 'http://hgdownload.soe.ucsc.edu/goldenPath/'+BUILD_ALT_ID+'/multiz100way/maf/',
                  '-O', args.maf_directory+'directory_list.html'])
        file_list = []
        with open(args.maf_directory+'directory_list.html') as directory_html:
            for html_line in directory_html:
                if html_line.strip().startswith('<a href="'):
                    file_list.append(html_line.strip().split('"')[1])

        for maf_file in file_list:
            sys.stderr.write('Downloading '+maf_file+'... ')
            call(['wget', 'http://hgdownload.soe.ucsc.edu/goldenPath/'+BUILD_ALT_ID+'/multiz100way/maf/'+maf_file,
                  '-O', args.maf_directory+'mafs/'+maf_file])
            sys.stderr.write('Done.')

    # ----------------------------------------------------------------------------------------------------
    elif args.create_exon_alignments:
        """Parse the MultiZ alignment files, by chromosome, and extract "gene-based" exon alignments"""
        # select the set of chromosomes to run on
        chroms = [('' if args.chromosome.startswith('chr') else 'chr') + args.chromosome]
        extract_exon_fasta_from_maf(chroms, args.maf_directory+'mafs/', args.aln_directory)

    # ----------------------------------------------------------------------------------------------------
    elif args.create_protein_alignments:
        """Stitch together the appropriate exon cDNA alignment files to generate 100-way protein alignments"""
        # select the set of chromosomes to run on
        chroms = [('' if args.chromosome.startswith('chr') else 'chr') + args.chromosome]
        create_protein_alignment(chroms, args.aln_directory)

    # ----------------------------------------------------------------------------------------------------
    elif args.compute_jsd:
        """Calculate JSD scores for the 100-way protein alignments in the corresponding directories"""
        # select the set of chromosomes to run on
        chroms = [('' if args.chromosome.startswith('chr') else 'chr') + args.chromosome]
        calculate_jsd_scores(chroms, args.aln_directory)

    # ----------------------------------------------------------------------------------------------------
    elif args.create_track_files:
        """Create two output files corresponding to the track intervals and track weights per protein"""
        track_file = args.outdir + '100way-jsdconservation_domsbyprot-' + GENOME_BUILD + '.txt.gz'
        weight_file = args.outdir + '100way-jsdconservation_domainweights-' + GENOME_BUILD + '.txt.gz'
        create_100way_jsd_files(track_file, weight_file, args.aln_directory)
