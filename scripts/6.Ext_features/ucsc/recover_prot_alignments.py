#!/usr/bin/python

"""
Map from UCSC Genome Browser identifiers to Ensembl identifiers in order to create FASTA-format
multiple sequence alignments
"""

import os
import sys
import gzip
from subprocess import call
from maf_to_fasta import calc_jsd_conservation

########################################################################################################

GENOME_BUILD = 'GRCh38'

PATH = '/home/snadimpa/workspace/cancertf/data/'
DBPATH = '/home/snadimpa/datadb/'
RESULTS_DIR = DBPATH+'ucscgb/hg38alignment/foranat/'
HUMAN_PROTFILE = DBPATH + 'ensembl/Homo_sapiens.'+GENOME_BUILD+'/Homo_sapiens.'+GENOME_BUILD+'.pep.all.verified.fa.gz'


########################################################################################################

def ucscgb_to_trans(mapping_file=DBPATH+'ucscgb/hg38alignment/ucscgb_to_ensembl.txt'):
  """
  :param mapping_file: full path to a tab-delimited file where column 0 is the name of a UCSC genome browser
                      protein ID, and column 11 is an Ensembl transcript ID, downloaded from
                      http://genome.ucsc.edu/cgi-bin/hgTables
  :return: dictionary of UCSC ID -> Ensembl transcript ID
  """

  ucscgb_to_ensembl = {}

  id_handle = gzip.open(mapping_file) if mapping_file.endswith('gz') else open(mapping_file)
  for id_line in id_handle:
    if id_line.startswith('#') or len(id_line[:-1].split('\t')) < 12:
      continue
    ucsc_id = id_line[:-1].split('\t')[0]
    align_id = id_line[:-1].split('\t')[11].split('.')[0]  # remove the version

    ucscgb_to_ensembl[ucsc_id] = align_id
  id_handle.close()

  return ucscgb_to_ensembl


########################################################################################################

def trans_to_protein(fasta_file=HUMAN_PROTFILE):
  """
  :param fasta_file: full path to a FASTA formatted file with chromosome, gene, protein, and transcript
                     information available in the header
  :return: dictionary of transcript ID -> (chromosome, gene, protein) IDs
  """

  trans_to_protloc = {}

  fasta_handle = gzip.open(fasta_file) if fasta_file.endswith('gz') else open(fasta_file)
  for fasta_line in fasta_handle:
    if fasta_line.startswith('>'):
      prot_id = fasta_line[fasta_line.find('prot:')+5:-1].split()[0].split('.')[0]
      gene_id = fasta_line[fasta_line.find('gene:')+5:-1].split()[0].split('.')[0]
      tran_id = fasta_line[fasta_line.find('transcript:')+11:-1].split()[0].split('.')[0]
      chromosome = fasta_line[fasta_line.find('chromosome:'+GENOME_BUILD+':'):].split()[0].split(':')[2]

      if tran_id not in trans_to_protloc:
        trans_to_protloc[tran_id] = set()
      trans_to_protloc[tran_id].add((chromosome, gene_id, prot_id))
  fasta_handle.close()

  return trans_to_protloc


########################################################################################################

def parse_alignment(alignment_file=DBPATH+'ucscgb/hg38alignment/knownGene.exonAA.fa.gz',
                    output_dir=RESULTS_DIR):
  """
  :param alignment_file: full path to a FASTA formatted multiple alignment file BY PROTEIN, downloaded from:
                         http://hgdownload.soe.ucsc.edu/goldenPath/hg38/multiz100way/alignments/knownGene.exonAA.fa.gz
  :param output_dir: full path to a directory to store chromosome -> gene ID -> protein alignments and JSD values
  :return: None, but print success messages
  """

  # retrieve the ID mappings:
  ucscid_to_tran = ucscgb_to_trans()
  tranid_to_protloc = trans_to_protein()

  # write out to files:
  current_files = None  # keep track of the file name we are writing to
  current_file_handles = None  # and the handle itself

  aln_handle = gzip.open(alignment_file) if alignment_file.endswith('gz') else open(alignment_file)
  for aln_line in aln_handle:
    if aln_line.startswith('>') and '_hg38_' in aln_line:

      # close the previous file that we had been writing to
      if current_files:
        for handle in current_file_handles:
          handle.close()

      # try to find the new file to write to:
      current_id = aln_line[1:-1].split('_')[0]
      if current_id in ucscid_to_tran and ucscid_to_tran[current_id] in tranid_to_protloc:
        current_files = []
        current_file_handles = []
        for chrom, gene_id, prot_id in tranid_to_protloc[ucscid_to_tran[current_id]]:
          if chrom not in os.listdir(output_dir):
            call(['mkdir', output_dir+chrom])
          if gene_id not in os.listdir(output_dir+chrom):
            call(['mkdir', output_dir+chrom+'/'+gene_id])
          current_files.append(output_dir+chrom+'/'+gene_id+'/'+prot_id+'.100way-alignment.fa.gz')

          current_file_handles.append(gzip.open(current_files[-1], 'a' if os.path.isfile(current_files[-1]) else 'w'))
      else:
        sys.stderr.write('No such ID: '+current_id+' '+ucscid_to_tran.get(current_id, '')+'\n')
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
  sys.stderr.write('Wrote to '+output_dir+'\n')


########################################################################################################

if __name__ == "__main__":

  if False:
    parse_alignment()

  if True:
    aln_suffix = '.100way-alignment.fa.gz'

    for chrom in os.listdir(RESULTS_DIR):
      if os.path.isdir(RESULTS_DIR+chrom):
        for gene in os.listdir(RESULTS_DIR+chrom):
          if os.path.isdir(RESULTS_DIR+chrom+'/'+gene):
            alns = [aln_file for aln_file in os.listdir(RESULTS_DIR+chrom+'/'+gene) if aln_file.endswith(aln_suffix)]

            for aln in alns:
              calc_jsd_conservation(RESULTS_DIR+chrom+'/'+gene+'/'+aln,
                                    RESULTS_DIR+chrom+'/'+gene+'/'+aln.replace(aln_suffix, '.jsd.txt'), '')

