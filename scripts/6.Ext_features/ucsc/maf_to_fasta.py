"""
Separate the 100-way vertebrate alignment(s) from UCEC's Genome Browser into separate directories
corresponding to the known protein-coding genes

Usage: maf_to_fasta.py <chromosome number>
"""

import os
import sys
import argparse
import gzip
import math
from subprocess import call
from bx.align import maf
from intervaltree import IntervalTree
from create_protein_trackfiles import mutational_likelihoods_overall

####################################################################################################
# CONSTANTS
####################################################################################################

GENOME_BUILD = 'GRCh38'  # 'GRCh37'
BUILD_ALT_ID = 'hg38'  # 'hg19'

if GENOME_BUILD == 'GRCh38':
  SPCS = ['GRCh38/hg38', 'CSAC 2.1.4/panTro4', 'gorGor3.1/gorGor3', 'WUGSC 2.0.2/ponAbe2', 'GGSC Nleu3.0/nomLeu3',
          'BGI CR_1.0/rheMac3', 'Macaca_fascicularis_5.0/macFas5', 'Baylor Panu_2.0/papAnu2',
          'Chlorocebus_sabeus 1.1/chlSab2', 'WUGSC 3.2/calJac3', 'Broad/saiBol1', 'Broad/otoGar3',
          'TupChi_1.0/tupChi1', 'Broad/speTri2', 'JacJac1.0/jacJac1', 'MicOch1.0/micOch1',
          'C_griseus_v1.0/criGri1', 'MesAur1.0/mesAur1', 'GRCm38/mm10', 'RGSC 6.0/rn6',
          'Broad HetGla_female_1.0/hetGla2', 'Broad/cavPor3', 'ChiLan1.0/chiLan1', 'OctDeg1.0/octDeg1',
          'Broad/oryCun2', 'OchPri3.0/ochPri3', 'SGSC Sscrofa10.2/susScr3', 'Vicugna_pacos-2.0.1/vicPac2',
          'CB1/camFer1',
          'Baylor Ttru_1.4/turTru2', 'Oorc_1.1/orcOrc1', 'PHO1.0/panHod1', 'UMD_3.1.1/bosTau8', 'ISGC Oar_v3.1/oviAri3',
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
          'Broad oreNil1.1/oreNil2', 'NeoBri1.0/neoBri1', 'AstBur1.0/hapBur1', 'MetZeb1.1/mayZeb1', 'PunNye1.0/punNye1',
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
          'Broad oreNil1.1/oreNil2', 'NeoBri1.0/neoBri1', 'AstBur1.0/hapBur1', 'MetZeb1.1/mayZeb1', 'PunNye1.0/punNye1',
          'NIG/UT MEDAKA1/oryLat2', 'Xiphophorus_maculatus-4.4.2/xipMac1', 'Broad/gasAcu1',
          'Genofisk GadMor_May2010/gadMor1', 'Zv9/danRer7', 'Astyanax_mexicanus-1.0.2/astMex1', 'LepOcu1/lepOcu1',
          'WUGSC 7.0/petMar2']
SPCS = [a.split('/')[-1] for a in SPCS]

CHROMOSOMES = map(lambda chr_val: 'chr' + str(chr_val), list(range(1, 23)) + ['X', 'Y'])

FULL_CHROM_IDS = {'chrKI270711.1': 'chr1_KI270711v1_random',
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

if 'Genomics/grid' in os.getcwd():
  DBPATH = '/Genomics/grid/users/snadimpa/data/'
  EXONPATH = '/Genomics/grid/users/snadimpa/data/ensembl/Homo_sapiens.' + GENOME_BUILD + '/'
  MAFPATH = '/Genomics/grid/users/snadimpa/data/ucscgb/' + BUILD_ALT_ID + 'alignment/'
  PROTFILE = '/Genomics/grid/users/snadimpa/data/ensembl/Homo_sapiens.' + GENOME_BUILD + '/Homo_sapiens.' + GENOME_BUILD + \
             '.pep.all.verified.fa.gz'
  runlocal = False

else:
  DBPATH = '/media/vineetb/t5-vineetb/dsprint/in/temp/'
  EXONPATH = '/home/snadimpa/datadb/ensembl/Homo_sapiens.' + GENOME_BUILD + '/'
  MAFPATH = '/home/snadimpa/datadb/ucscgb/' + BUILD_ALT_ID + 'alignment/'
  PROTFILE = '/home/snadimpa/datadb/ensembl/Homo_sapiens.' + GENOME_BUILD + '/Homo_sapiens.' + GENOME_BUILD + \
             '.pep.all.verified.fa.gz'
  runlocal = True


####################################################################################################

def get_exon_locs(protfile=PROTFILE, chrom=CHROMOSOMES):
  """
  :param gene_path: Path to directory of Ensembl genes
  :return: an IntervalTree of gene locations on the appropriate chromosome
  """

  genelocs = {chromID: IntervalTree() for chromID in chrom}

  if not os.path.isfile(protfile):
    sys.stderr.write('Could not open ' + protfile + '\n')
    return genelocs

  x = gzip.open(protfile) if protfile.endswith('gz') else open(protfile)

  locstart = 'chromosome:' + GENOME_BUILD + ':'
  geneintervals = {}

  for l in x:
    if l.startswith('>'):
      loc = l[l.find(locstart) + len(locstart):].split()[0]
      chromID = loc.split(':')[0]
      chromID = ('' if chromID.startswith('chr') else 'chr') + chromID

      if chromID not in chrom: continue

      exons = loc.split(':')[1].replace('join(', '').replace('complement(', '').replace(')', '').split(',')
      exons = [tuple(sorted([int(a.split('..')[0]), int(a.split('..')[1])])) for a in exons]

      geneID = l[l.find('gene:') + 5:].split()[0]  # e.g., ENSG00000000457.9

      if chromID not in geneintervals:
        geneintervals[chromID] = {}
      if geneID not in geneintervals[chromID]:
        geneintervals[chromID][geneID] = set()

      # Add each POSITIVE exon SEPARATELY!
      for exon_start, exon_end in exons:  # these are 1-indexed, but we need 0-indexed...
        if exon_start > -1 and exon_end > -1:
          geneintervals[chromID][geneID].add((exon_start, exon_end))
  x.close()

  # Store all the information as interval trees:
  for chromID in geneintervals.keys():
    for geneID in geneintervals[chromID].keys():
      for exon_start, exon_end in geneintervals[chromID][geneID]:
        genelocs[chromID].addi(exon_start, exon_end + (1 if exon_end == exon_start else 0),
                               geneID + '_' + str(exon_start) + '-' + str(exon_end))

  return genelocs


####################################################################################################

def concatenate_maf_to_fasta(maffile, fastafile, startpos, endpos):
  """
  :param maffile: full path to a MultiZWay alignment file
  :param fastafile: full path to resulting fasta multiple alignment file
  :return: None
  """

  species = SPCS  # all the species that should be represented in file
  fill = ""  # how to separate blocks?
  wrap = -1  # maximum line length in alignment (or no maximum)

  OUT = gzip.open(fastafile, 'w') if fastafile.endswith('gz') else open(fastafile, 'w')

  texts = {}  # store the complete sequence for each species
  for s in species:
    texts[s] = []

  # Read in the maf file !
  maf_reader = maf.Reader(gzip.open(maffile) if maffile.endswith('gz') else open(maffile))

  for m in maf_reader:
    for s in species:
      c = m.get_component_by_src_start(s)
      if c:
        texts[s].append(c.text)
      else:
        texts[s].append("-" * m.text_size)

  # Finally, print out the concatenated sequence by species in order
  for s in species:
    OUT.write('>' + s + (' ' + str(startpos + 1) + ':' + str(endpos + 1) if s == BUILD_ALT_ID else '') + '\n')
    seqmatch = fill.join(texts[s])
    if wrap <= 0:
      OUT.write(seqmatch + '\n')

    else:
      p = 0  # starting index of entire alignment for each line
      while p < len(seqmatch):
        OUT.write(seqmatch[p:min(p + wrap, len(seqmatch))] + '\n')
        p += wrap

  OUT.close()
  sys.stderr.write('Finished!\n')


####################################################################################################

def translate_rna(rnaSeq):
  """
  translate the cDNA sequence to protein using GENCODE; unknown codons
  or accidental stop codons are translated as 'X'
  """

  protSeq = []
  for i in range(0, len(rnaSeq), 3):
    codon = rnaSeq[i:i + 3]
    if len(codon) == 3:
      if codon not in GENCODE and '-' in codon:
        protSeq.append('-')
      elif codon not in GENCODE and 'N' in codon:
        protSeq.append('X')
      elif GENCODE[codon] == "_" and i < len(rnaSeq) - 5:  # stop codon somewhere in the middle..
        protSeq.append('X')
      else:
        protSeq.append(GENCODE[codon])
  return ''.join(protSeq)


####################################################################################################

def stitchexons(exonfile, exon_alignment_files, aln_outfile, error_check=False):
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
  negcomplement = False  # are we looking at the complementary strand for the ENTIRE PROTEIN?

  exon_handle = gzip.open(exonfile) if exonfile.endswith('gz') else open(exonfile)
  for l in exon_handle:
    if l.startswith('>'):
      if 'complement' in l:
        negcomplement = True
      continue

    complement = False  # is this particular exon already reversed?
    exonloc, exonseq = l[:-1].split('\t')

    exonstart, exonend = map(int, exonloc.split(':'))

    # swap the start and end if necessary
    if exonend < exonstart or (exonend < 0 and exonstart < 0 and exonend == exonstart and negcomplement):
      complement = True
      tmp = exonstart
      exonstart = exonend
      exonend = tmp

    # if we can't find the corresponding exon alignment, fluff this exon with gaps, and move on
    if (str(exonstart), str(exonend)) not in exon_alignment_files:
      sys.stderr.write('No alignment file for ' + exonfile + ': ' + str(exonstart) + '-' + str(exonend) + '\n')
      cdna_per_spc[BUILD_ALT_ID].append(''.join([COMPLEMENT[i] for i in list(exonseq)]) if complement else exonseq)
      for spc in SPCS[1:]:
        cdna_per_spc[spc].append('-' * len(exonseq))
      continue

    # otherwise, keep track of the exon information !
    with gzip.open(exon_alignment_files[(str(exonstart), str(exonend))]) as exon_aln_handle:
      human_block_location = exon_aln_handle.next()[:-1].split()[1]
      blockstart, blockend = map(int, human_block_location.split(':'))  # this is ALWAYS in the correct order

      humanseq = exon_aln_handle.next()[:-1]

      # these indices will correspond to the indices into the sequence with gaps (i.e. not continuous)
      nongap_indices = [(nt, i) for (nt, i) in [(nt, i) for i, nt in enumerate(list(humanseq))] if nt != '-']

      if error_check:
        print('------------------')
        print('exonstart = ' + str(exonstart))
        print('blockstart = ' + str(blockstart))
        print('exonend = ' + str(exonend))
        print('blockend = ' + str(blockend))

      # if the block is SMALLER than our desired exon, pad out the block accordingly:
      fluffleft = 0 if exonstart >= blockstart else blockstart - exonstart + 1
      fluffright = 0 if exonend <= blockend else exonend - blockend + 1

      if complement:
        tmp = fluffleft
        fluffleft = fluffright
        fluffright = tmp

      # shave off excess nucleotides from front and back:
      start_exon_index = 0 if exonstart <= blockstart else exonstart - blockstart  # shave off excess in front
      end_exon_index = len(nongap_indices) if exonend >= blockend else len(nongap_indices) - (blockend - exonend - 1)

      if error_check:
        print('fluffleft: ' + str(fluffleft))
        print('fluffright: ' + str(fluffright))
        print('start_block_index: ' + str(start_exon_index))
        print('end_block_index: ' + str(end_exon_index))

      thisexon = nongap_indices[start_exon_index:end_exon_index]
      indices = [index for nuc, index in thisexon]
      if complement:
        indices = indices[::-1]

      fixed = True
      if ''.join([humanseq[i] for i in indices]).upper() != exonseq.upper()[fluffleft:len(exonseq) - fluffright]:

        if error_check:
          print('us-----: ' + ''.join([humanseq[i] for i in indices]).upper())
          print('correct: ' + exonseq.upper()[fluffleft:len(exonseq) - fluffright])

        fixed = False

        newexons = {'remove left': nongap_indices[start_exon_index + 1:end_exon_index],
                    'remove right': nongap_indices[start_exon_index:end_exon_index - 1],
                    'add left': nongap_indices[start_exon_index - (1 if start_exon_index > 0 else 0):end_exon_index],
                    'add right': nongap_indices[
                                 start_exon_index:end_exon_index + (1 if end_exon_index < len(nongap_indices) else 0)],
                    'shift left': nongap_indices[
                                  start_exon_index - (1 if start_exon_index > 0 else 0):end_exon_index - 1],
                    'shift right': nongap_indices[start_exon_index + 1:end_exon_index + (
                      1 if end_exon_index < len(nongap_indices) else 0)]}

        for fixtype, newexon in newexons.items():
          fixed_indices = [index for nuc, index in newexon]
          if complement:
            fixed_indices = fixed_indices[::-1]

          # did we try and fail to add left?
          orig_fluff_left = fluffleft
          orig_fluff_right = fluffright
          if fixtype in {'add left', 'shift left'} and start_exon_index <= 0:
            if complement:
              fluffright += 1
            else:
              fluffleft += 1
          if fixtype in {'add right', 'shift right'} and end_exon_index >= len(nongap_indices):
            if complement:
              fluffleft += 1
            else:
              fluffright += 1

          if error_check:
            print(fixtype.split()[0][0] + fixtype.split()[1][0] + '-----: ' + \
                  ''.join([humanseq[i] for i in fixed_indices]).upper())

          if ''.join([humanseq[i] for i in fixed_indices]).upper() == exonseq.upper()[
                                                                      fluffleft:len(exonseq) - fluffright]:
            indices = fixed_indices
            fixed = True
            cfixtype = fixtype if not complement else (fixtype.replace('left', 'right') if 'left' in fixtype else
                                                       fixtype.replace('right', 'left'))
            sys.stderr.write('Fixed by ' + cfixtype + '\n')
            break

          fluffleft = orig_fluff_left
          fluffright = orig_fluff_right

      if not fixed:
        sys.stderr.write('Exon ' + exonfile + ': ' + str(exonstart) + '-' + str(exonend) + ' did not match:\n' +
                         ''.join([humanseq[i] for i in indices]).upper() + '\n' +
                         exonseq.upper() + '\n')

        for fixtype, newexon in newexons.items():
          cfixtype = fixtype if not complement else (fixtype.replace('left', 'right') if 'left' in fixtype else
                                                     fixtype.replace('right', 'left'))
          sys.stderr.write(cfixtype + '\n' +
                           ''.join([humanseq[i] for i in [a[1] for a in newexon][::-1 if complement else 1]]).upper() +
                           '\n#\n')

        cdna_per_spc[BUILD_ALT_ID].append(''.join([COMPLEMENT[i] for i in list(exonseq)]) if complement else exonseq)
        for spc in SPCS[1:]:
          cdna_per_spc[spc].append('-' * len(exonseq))
        continue

      cdna_per_spc[BUILD_ALT_ID].append(('-' * fluffleft) + ''.join([COMPLEMENT[humanseq[i]] for i in indices])
                                        if complement else ''.join([humanseq[i] for i in indices]) + ('-' * fluffright))

      for aln_line in exon_aln_handle:
        if aln_line.startswith('>'):
          currspc = aln_line[1:-1]
          currseq = exon_aln_handle.next()[:-1]
          cdna_per_spc[currspc].append(('-' * fluffleft) + ''.join(COMPLEMENT[currseq[i]] for i in indices)
                                       if complement else ''.join(currseq[i] for i in indices) + ('-' * fluffright))
  exon_handle.close()

  out_handle = gzip.open(aln_outfile, 'w') if aln_outfile.endswith('gz') else open(aln_outfile, 'w')
  for spc in SPCS:
    # out_handle.write('>' + spc + '\n' + ''.join(cdna_per_spc[spc]).upper() + '\n')
    out_handle.write('>' + spc + '\n' + translate_rna(''.join(cdna_per_spc[spc]).upper()) + '\n')
    out_handle.close()
  sys.stderr.write('Wrote multiple protein alignment to ' + aln_outfile + '\n')


####################################################################################################

def fix_cdna_prots(chrom, gene, prot):
  """
  Given the exon file, recreate the cDNA and protein files.... (ALWAYS chop off last character of protein!)
  """

  exonfile = EXONPATH + chrom + '/' + gene + '/' + prot + '.exons.txt'
  cdnafile = open(EXONPATH + chrom + '/' + gene + '/' + prot + '.cdna.fa', 'w')
  protfile = open(EXONPATH + chrom + '/' + gene + '/' + prot + '.prot.fa', 'w')

  header = ''
  rnaseq = []

  with open(exonfile) as x:
    header = x.next()
    complement = 'complement' in header

    for l in x:
      exon = l[:-1].split()[1]
      rnaseq.append(''.join([COMPLEMENT[i] for i in list(exon)]) if complement else exon)

  cdnafile.write(header + ''.join(rnaseq) + '\n')
  cdnafile.close()

  protfile.write(header + translate_rna(''.join(rnaseq).upper())[:-1] + '\n')
  protfile.close()


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
  x = gzip.open(alnfile) if alnfile.endswith('gz') else open(alnfile)
  for l in x:
    if l.startswith('>'):
      spc = l[1:-1].split()[0].split('_')[1]  # species name (e.g., hg38...)
      if spc not in seqs:
        seqs[spc] = []
      if mainspc == '':
        mainspc = spc
      seqs[spc].append(x.next().strip())  # full protein sequence (NOTE: ends with stop codon)
  x.close()

  for spc in seqs.keys():
    seqs[spc] = ''.join(seqs[spc])
    if seqs[spc].endswith('Z'):
      seqs[spc] = seqs[spc][:-1]  # remove "stop" codon

  # NOW, calculate the column-based JSD...
  aa = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S',
        'T', 'W', 'Y', 'V']  # 20 amino acids, not including gaps...

  # Output file containing the JSD conservation values for each position
  OUT = open(outfile, 'w')
  OUT.write(header)

  for i in xrange(len(seqs[mainspc])):
    caa = seqs[mainspc][i]  # residue value at this position in the human sequence
    ccol = [cseq[i] for cseq in seqs.values() if i < len(cseq) and cseq[i] in aa]  # full ungapped column
    cval = len(ccol) / 100. * JSdiv_column(ccol)  # overall conservation score
    OUT.write('\t'.join([str(i), str(caa), str(cval)]) + '\n')
  OUT.close()
  sys.stderr.write('Wrote ' + outfile + '\n')


####################################################################################################

def check_alignments(chrom, gene, prot):
  """
  :param chrom: chromosome ID
  :param gene: gene ID
  :param prot: protein ID
  :return: whether the alignment sequence matches (well enough) the expected protein sequence
  """

  alignment_file = EXONPATH + chrom + '/' + gene + '/' + prot + '.100way-alignment.fa.gz'
  protein_file = EXONPATH + chrom + '/' + gene + '/' + prot + '.prot.fa'

  with gzip.open(alignment_file) as aln_file:
    aln_file.next()
    aln_seq = aln_file.next().strip()

  with open(protein_file) as prot_file:
    prot_file.next()
    prot_seq = prot_file.next().strip()

  distance_allowance = 0.05 * len(prot_seq)
  hamming_distance = sum(
    [1 if prot_seq[i] != (aln_seq[i] if i < len(aln_seq) else '-') else 0 for i in xrange(len(prot_seq))])

  if hamming_distance > distance_allowance:
    print(EXONPATH + chrom + '/' + gene + '/' + prot)
    return False
  else:
    return True


####################################################################################################

def create_exac_jsd_files(domainfile=DBPATH + 'exac/exac-jsdconservation_domsbyprot-' + GENOME_BUILD + '.txt.gz',
                          weightfile=DBPATH + 'exac/exac-jsdconservation_domainweights-' + GENOME_BUILD + '.txt.gz'):
  """
  :param domainfile: full path to a tab-delimited file containing the LOCATIONS of domains (match state to
                     0-indexed protein position
  :param weightfile: full path to a tab-delimited file containing the BINDING POTENTIAL SCORES of each match
                     state in a domain
  :return: None, but print success message
  """

  weight_handle = gzip.open(weightfile, 'w') if weightfile.endswith('gz') else open(weightfile, 'w')
  weight_handle.write('\n'.join(
    ['# Continuous weights across entire protein sequences corresponding to their ExAC minor allele frequency',
     '#   ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz',
     '#   Information downloaded on January 25, 2017',
     '# MAF values are binned as follows: ',
     '\n'.join(['#   1 if 0 < AF <= 0.001',
                '#   0.75 if 0.001 < AF <= .005',
                '#   0.5 if 0.005 < AF <= 0.01',
                '#   0.25 if 0.01 < AF <= 0.1',
                '#   0 if 0.1 < AF']),
     '\t'.join(['#EnsemblProtID_ExAC_JSD', 'LigandTypeAlways-ExAC_CONS_', 'MatchState',
                'ConservationScore', 'NumInstancesAlways10',
                '1-Indexed-ProteinPositionAndResidue']) + '\n']))

  dom_handle = gzip.open(domainfile, 'w') if domainfile.endswith('gz') else open(domainfile, 'w')
  dom_handle.write('\n'.join(['# All proteins for which at least one position contained a natural variant in ExAC',
                              '#   ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/subsets/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz',
                              '#   Information downloaded on January 25, 2017',
                              '\t'.join(['#EnsemblProtID', 'DomainID', 'matchstate:AA-0-index:AA-value\n'])]))

  # Now, process the input file..
  # protein_id -> (protein_length, ((position, AF), (position, AF), (position, AF),...))
  all_mutation_rates = get_original_background_mutation(DBPATH + 'tcga/exac_allelefreqs-' + GENOME_BUILD + '.txt.gz')

  # Get the protein sequences for proteins that we are able to model:
  id_to_seq = {}
  seq_handle = gzip.open(PROTFILE) if PROTFILE.endswith('gz') else open(PROTFILE)
  for seq_line in seq_handle:
    if seq_line.startswith('>') and seq_line[1:-1].split()[0] in all_mutation_rates:
      prot_id = seq_line[1:-1].split()[0]
      if prot_id in all_mutation_rates:
        id_to_seq[prot_id] = seq_handle.next().strip()
  seq_handle.close()

  # Finally, write out the information we need:
  for prot_id, seq in sorted(id_to_seq.items()):
    position_to_weight = {}

    # pos is the 0-indexed protein position (int)
    for pos, af in all_mutation_rates[prot_id][1]:

      if af == 0:  # we will only model sites with some sort of allele frequency
        continue

      position_to_weight[pos] = (1.0 if af <= .001 else
                                 (0.75 if af <= .005 else
                                  (0.5 if af <= .01 else
                                   (0.25 if af <= .1 else 0.))))
    match_states = []
    mstate_index = 0
    for i, aa in enumerate(list(seq)):
      if i in position_to_weight:
        match_states.append(':'.join([str(mstate_index + 1), str(i), aa]))  # match state, sequence index, amino acid
        weight_handle.write('\t'.join([prot_id + '_ExAC_JSD', 'ExAC_CONS_', str(mstate_index + 1),
                                       str(position_to_weight[i]), '10', str(i + 1) + ':' + aa]) + '\n')
        mstate_index += 1
    dom_handle.write('\t'.join([prot_id, prot_id + '_ExAC_JSD', ','.join(match_states)]) + '\n')

  dom_handle.close()
  sys.stderr.write('Wrote to ' + domainfile + '\n')

  weight_handle.close()
  sys.stderr.write('Wrote to ' + weightfile + '\n')


####################################################################################################

def create_100way_jsd_files(domainfile=MAFPATH + '100way-jsdconservation_domsbyprot-' + GENOME_BUILD + '.txt.gz',
                            weightfile=MAFPATH + '100way-jsdconservation_domainweights-' + GENOME_BUILD + '.txt.gz'):
  """
  Create two types of files: one (domsbyprot) and the second (biolipsummary).
  """

  weightfile_handle = gzip.open(weightfile, 'w') if weightfile.endswith('gz') else open(weightfile, 'w')
  weightfile_handle.write('\n'.join(
    ['# Continuous conservation weights across entire protein sequences, calculated using UCSC Genome Browser\'s',
     '#   100-vertebrate multiple alignment (http://hgdownload.soe.ucsc.edu/goldenPath/' + BUILD_ALT_ID + '/multiz100way/), downloaded June 28, 2017.',
     '# Nucleotide alignments for all protein isoforms were extracted and translated to create per-protein multiple alignments.',
     '# Conservation scores per column were calculated as the Jensen-Shannon divergence of the non-gapped residues',
     '#   against a Blosum 62 background (https://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt), multiplied by',
     '#   the fraction of non-gapped residues in the column.',
     '\t'.join(['#EnsemblProtID_JSD', 'LigandTypeAlwaysCONS_', '1-Indexed-ProteinPosition', 'ConservationScore',
                'NumInstancesAlways10', '1-Indexed-ProteinPositionAndResidue']) + '\n']))

  domainfile_handle = gzip.open(domainfile, 'w') if domainfile.endswith('gz') else open(domainfile, 'w')
  domainfile_handle.write(
    '\n'.join(['# All proteins for which conservation scores could be determined, using the UCSC Genome Browser\'s',
               '#   100-vertebrate multiple alignment (http://hgdownload.soe.ucsc.edu/goldenPath/' + BUILD_ALT_ID +
               '/multiz100way/), downloaded June 28, 2017.',
               '\t'.join(['#EnsemblProtID', 'DomainID', 'matchstate:AA-0-index:AA-value\n'])]))

  for chrom in sorted(os.listdir(EXONPATH)):
    if not os.path.isdir(EXONPATH + chrom): continue

    for geneID in sorted(os.listdir(EXONPATH + chrom)):
      if not os.path.isdir(EXONPATH + chrom + '/' + geneID): continue

      for protID in sorted(
          [a.replace('.jsd.txt', '') for a in os.listdir(EXONPATH + chrom + '/' + geneID) if a.endswith('.jsd.txt')]):
        conservation_score_available = False
        fullprotein = []

        for l in open(EXONPATH + chrom + '/' + geneID + '/' + protID + '.jsd.txt'):
          if l.startswith('>'): continue
          protpos, residue, cons = l.strip().split('\t')[:3]
          if len(protpos) < 1 or len(residue) < 1 or len(cons) < 1: continue

          cons = max(0, float(cons))
          fullprotein.append(str(int(protpos) + 1) + ':' + protpos + ':' + residue)

          if cons > 0:
            conservation_score_available = True
            weightfile_handle.write('\t'.join([protID + '_100way_JSD', '100WAY_CONS_', str(int(protpos) + 1), str(cons),
                                        '10', str(int(protpos) + 1) + ':' + residue]) + '\n')
        if conservation_score_available:
          # last "amino acid" is SOMETIMES a stop codon, so remove it to prevent failures
          if fullprotein[-1].split(':')[2] == '_': fullprotein = fullprotein[:-1]
          domainfile_handle.write('\t'.join([protID, protID + '_100way_JSD', ','.join(fullprotein)]) + '\n')

  domainfile_handle.close()
  weightfile_handle.close()


####################################################################################################
# MAIN
####################################################################################################

if __name__ == "__main__":

  parser = argparse.ArgumentParser(description='Convert multiple alignment into separate gene-based files.')

  parser.add_argument('--chromosome', type=str,
                      help='Chromosome number in ' + GENOME_BUILD + ' (i.e., 1-22, X, Y).',
                      default='chr1')
  args = parser.parse_args()

  if not runlocal:
    chroms = [('' if args.chromosome.startswith('chr') else 'chr') + args.chromosome]
  else:
    chroms = CHROMOSOMES

  # ----------------------------------------------------------------------------------------------------
  if True:
    """
    Create "alignments" using the ExAC minor allele frequencies in order to compute column-wise JSDs
    """
    create_exac_jsd_files()
    sys.exit(1)


  # ----------------------------------------------------------------------------------------------------
  if False:
    """
    Outdated code to generate list of domain matches (this has been migrated to runhmmer.py, as well as
    code to use the already-processed 100-way alignments and determine the JSDs per protein position
    """
    create_100way_jsd_files()
    sys.exit(1)


  # ----------------------------------------------------------------------------------------------------
  if False:
    """
    Stitch together the appropriate exon cDNA alignment files to generate 100-way protein alignments
    """
    total_failed = 0
    total_processed = 0
    for chrom_id in chroms:
      if not os.path.isdir(EXONPATH + chrom_id.replace('chr', '')):
        continue

      geneIDs = os.listdir(EXONPATH + chrom_id.replace('chr', '') + '/')

      for geneID in geneIDs:
        if not os.path.isdir(EXONPATH + chrom_id.replace('chr', '') + '/' + geneID + '/'): continue

        for protID in [a.replace('.exons.txt', '') for a in
                       os.listdir(EXONPATH + chrom_id.replace('chr', '') + '/' + geneID)
                       if a.endswith('.exons.txt')]:

          # full path to the exons file (sorted exons)
          exonfile = EXONPATH + chrom_id.replace('chr', '') + '/' + geneID + '/' + protID + '.exons.txt'
          if not os.path.isfile(exonfile):
            continue

          header = ''
          with open(exonfile) as x:
            header = x.next()

          # set of all exon alignments (we'll see which ones we need for this particular protein)
          alignments = {(a.split('_')[1].split('-')[0],
                         a.split('-')[1].split('.')[0]): EXONPATH + chrom_id.replace('chr', '') + '/' + geneID + '/' + a
                        for a in os.listdir(EXONPATH + chrom_id.replace('chr', '') + '/' + geneID + '/') if
                        a.startswith(geneID) and a.endswith('.aln.fa.gz')}

          # output alignment file
          alignment_file = EXONPATH + chrom_id.replace('chr',
                                                      '') + '/' + geneID + '/' + protID + '.100way-alignment.fa.gz'

          # output JSD file:
          jsd_file = EXONPATH + chrom_id.replace('chr', '') + '/' + geneID + '/' + protID + '.jsd.txt'

          # fix_cdna_prots(chrom_id.replace('chr',''), geneID, protID)
          # stitchexons(exonfile, alignments, alignment_file)
          calc_jsd_conservation(alignment_file, jsd_file, header)
          # passing = check_alignments(chrom_id.replace('chr',''), geneID, protID)
          # total_processed += 1
          # if not passing:
          #  total_failed += 1

    print(str(total_failed) + '/' + str(total_processed))
    sys.exit(1)


  # ----------------------------------------------------------------------------------------------------
  if False:
    """
    Copy alignment files locally (to move off of cluster) so that we don't use up too much space...
    """

    for chrom_id in chroms:
      sys.stderr.write('Getting exon locations for genes on chromosome(s) ' + ','.join([chrom_id]) + '\n')
      gene_locs = get_exon_locs(PROTFILE, [chrom_id])  # All gene locations on this chromosome
      sys.stderr.write('Finished!\n')

      full_chrom_name = FULL_CHROM_IDS.get(chrom_id, chrom_id)

      maf_file = MAFPATH + chrom_id + '.maf.gz'
      gene_path = EXONPATH + chrom_id.replace('chr', '') + '/'

      if not os.path.isfile(maf_file):
        sys.stderr.write('Could not open ' + maf_file + '\n')
        if not runlocal:
          call(['scp', 'snadimpa@shilpa.princeton.edu:~/datadb/ucscgb/hg38alignment/' + chrom_id + '.maf.gz', maf_file])
        else:
          continue

      if not os.path.isdir(gene_path):
        sys.stderr.write('No such directory ' + gene_path + '\n')
        os.system('mkdir ' + gene_path)

      newfiles = {}  # The new smaller output .maf files and the regions they cover

      current_block = []  # Temporary list to keep track of each alignment chunk

      x = gzip.open(maf_file) if maf_file.endswith('gz') else open(maf_file)

      maf_header = []  # Each .maf file must start with ##maf ...
      for l in x:
        if not l.startswith('a score'):
          maf_header.append(l)
        else:
          current_block.append(l)
          break

      sys.stderr.write('Processed ' + maf_file + ' header...\n')

      for l in x:

        # Write out the current block:
        if l.startswith('\n') and len(current_block) > 0:
          current_block.append(l)
          startloc, seqlen = map(int,
                                 [a for a in current_block if BUILD_ALT_ID + '.' + full_chrom_name in a][0].split()[
                                 2:4])

          # Remember that the .maf files are 0-indexed whereas the exon locations are 1-indexed.

          for start, end, exonID in gene_locs[chrom_id][startloc + 1:startloc + seqlen + 1]:

            geneID = exonID.split('_')[0]

            genefile = gene_path + geneID + '/' + exonID + '.alignment.maf'

            if genefile not in newfiles:
              OUT = open(genefile, 'w')
              map(OUT.write, maf_header)
              newfiles[genefile] = [startloc, startloc + seqlen]
            else:
              OUT = open(genefile, 'a')
              newfiles[genefile][1] = startloc + seqlen

            map(OUT.write, current_block)
            OUT.close()
          current_block = []

        else:
          current_block.append(l)
      x.close()

      # Finally, convert to a fasta file !
      for exon_maf_file in sorted(list(newfiles.keys())):
        sys.stderr.write('Concatenating ' + exon_maf_file + ' to fasta...\n')
        concatenate_maf_to_fasta(exon_maf_file,
                                 exon_maf_file.replace('.alignment.maf', '.aln.fa.gz'),
                                 newfiles[exon_maf_file][0],
                                 newfiles[exon_maf_file][1])
        os.system('rm ' + exon_maf_file)

      if not runlocal:
        call(['rm', maf_file])
