#!/usr/bin/python

"""
Extract all Ensembl protein information (for all isoforms) and cross-check against the expected CDS and DNA
files; fix isoforms where necessary (e.g., some exons are off by 1-2 nucleotides).

Written by Dario Ghersi (dghersi@princeton.edu)
Edited by Pawel Przytycki (pprzytyc@princeton.edu)
Contact Shilpa N. Kobren (snadimpa@princeton.edu) with questions
"""

import os
import re
import sys
import collections
import gzip
import argparse
from datetime import date
from subprocess import call
from config import data_path, GENOME_BUILD as BUILD

####################################################################################################
# PARAMETERS (to update)
####################################################################################################

# Set this parameter to limit to chromosomes 1-22, X, Y to exclude scaffolds (speed things up)
VALID_CHROMOSOMES = None if True else map(str, range(1, 23)) + ['X', 'Y']

MAX_MISS_RATE = .1  # maximum mismatch between translated protein and expected protein for "match"

####################################################################################################
# CONSTANTS
####################################################################################################

TODAY = str(date.today().year).zfill(4) + \
        str(date.today().month).zfill(2) + \
        str(date.today().day).zfill(2)

HUMAN_PROTFILE = data_path + 'ensembl/Homo_sapiens.' + BUILD + '/Homo_sapiens.' + BUILD + '.pep.all.fa.gz'

PROGRESS_FILE = data_path + 'ensembl/Homo_sapiens.' + BUILD + '/exons/proteingeneration-' + TODAY + '.log'

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N",
              "a": "t", "t": "a", "c": "g", "g": "c", "n": "n"}

MAX_SEQ_LINE = 70

NONDECIMAL = re.compile(r'[^\d]+')


####################################################################################################
# Preliminaries: Update original FASTA file with additional information; inflate DNA files
####################################################################################################

def ensembl_to_altids(hgncfile=data_path+'ensembl/Homo_sapiens.'+BUILD+'/Homo_sapiens.'+BUILD+'.toHGNC.tsv',
                      refseqfile=data_path+'ensembl/Homo_sapiens.'+BUILD+'/Homo_sapiens.'+BUILD+'.toRefSeq.tsv',
                      entrezfile=data_path+'ensembl/Homo_sapiens.'+BUILD+'/Homo_sapiens.'+BUILD+'.toEntrez.tsv',
                      locfile=data_path+'ensembl/Homo_sapiens.'+BUILD+'/Homo_sapiens.'+BUILD+'.genelocs.tsv'):
    """
    :param hgncfile: full path to a tab-delimited file with columns Ensembl Gene ID, Ensembl Protein ID, HGNC ID,
                     Hugo Symbol, and Gene Name
    :param refseqfile: full path to a tab-delimited file with columns Ensembl Protein ID, RefSeq Protein ID,
                       RefSeq Predicted Protein ID
    :param entrezfile: full path to a tab-delimited file with columns Ensembl Protein ID, Entrez Gene ID
    :param locfile: full path to a tab-delimited file with columns Ensembl Gene ID, Ensembl Transcript ID,
                    Ensembl Protein ID, Chromosome ID, Strand, Exon Rank in Transcript, Genomic Exon Start,
                    Genomic Exon End
    :return: dictionaries of Ensembl Protein IDs -> alternate identifiers
    """

    # (1) Process HGNC IDs and Hugo Symbols:
    to_hgnc_id, to_hugo_symbol = {}, {}
    hgnc_inhandle = gzip.open(hgncfile) if hgncfile.endswith('gz') else open(hgncfile)
    hgnc_inhandle.next()  # skip the header
    for altid_line in hgnc_inhandle:
        ensembl_gene, ensembl_prot, hgnc_id, hugo_symbol, gene_name = altid_line[:-1].split('\t')[:5]

        # Set the HGNC ID if it is a valid ID
        if hgnc_id != '':
            if ensembl_prot != '':
                to_hgnc_id[ensembl_prot] = hgnc_id.replace('HGNC:', '')
            if ensembl_gene != '':
                to_hgnc_id[ensembl_gene] = hgnc_id.replace('HGNC:', '')

        # Set the Hugo Symbol if it's valid
        if hugo_symbol != '':
            if ensembl_prot != '':
                to_hugo_symbol[ensembl_prot] = hugo_symbol
            if ensembl_gene != '':
                to_hugo_symbol[ensembl_gene] = hugo_symbol
        elif gene_name != '':  # only IF hugo symbol is blank!
            if ensembl_prot != '':
                to_hugo_symbol[ensembl_prot] = gene_name
            if ensembl_gene != '':
                to_hugo_symbol[ensembl_gene] = gene_name
    hgnc_inhandle.close()

    # (2) Process RefSeq IDs
    to_refseq = {}
    refseq_inhandle = gzip.open(refseqfile) if refseqfile.endswith('gz') else open(refseqfile)
    refseq_inhandle.next()
    for altid_line in refseq_inhandle:
        if not altid_line.strip().startswith('ENSP') or len(altid_line.strip().split()) < 2:
            continue
        # this automatically excludes the PREDICTED ID, if there is one
        ensembl_prot, refseq_id = altid_line.strip().split()[:2]
        to_refseq[ensembl_prot] = refseq_id
    refseq_inhandle.close()

    # (3) Process Entrez IDs
    to_entrez = {}
    entrez_inhandle = gzip.open(entrezfile) if entrezfile.endswith('gz') else open(entrezfile)
    entrez_inhandle.next()
    for altid_line in entrez_inhandle:
        ensembl_prot, entrez_id = altid_line[:-1].split('\t')[:2]
        if ensembl_prot != '' and entrez_id != '':
            to_entrez[ensembl_prot] = entrez_id
    entrez_inhandle.close()

    # (4) Process Ensembl Gene and Transcript IDs *AND EXON LOCATIONS*
    to_gene, to_transcript, to_location = {}, {}, {}

    exons = {}
    ensembl_inhandle = gzip.open(locfile) if locfile.endswith('gz') else open(locfile)
    ensembl_inhandle.next()
    for altid_line in ensembl_inhandle:
        if len(altid_line.strip().split('\t')) < 8:
            continue
        ensembl_gene, ensembl_trans, ensembl_prot, chrom, strand, rank, start, end = altid_line[:-1].split('\t')[:8]

        if ensembl_prot not in exons:
            exons[ensembl_prot] = set()

        exons[ensembl_prot].add((int(rank), int(start), int(end), chrom, strand))
        to_gene[ensembl_prot] = ensembl_gene
        to_transcript[ensembl_prot] = ensembl_trans
    ensembl_inhandle.close()

    # standard location format is comma-delimited, ordered exons with genomic start/end locations specified
    for ensembl_prot, allx in sorted(exons.items()):
        loc = ','.join([str(exon_start) + '..' + str(exon_end) for _, exon_start, exon_end, _, _ in sorted(list(allx))])
        if ',' in loc:
            loc = 'join(' + loc + ')'
        if '-' in list(allx)[0][4]:  # strand is negative?
            loc = 'complement(' + loc + ')'
        to_location[ensembl_prot] = BUILD + ':' + list(allx)[0][3] + ':' + loc

    return to_hgnc_id, to_hugo_symbol, to_refseq, to_entrez, to_gene, to_transcript, to_location


####################################################################################################

def update_fasta_genelocs_altids(fastafile=HUMAN_PROTFILE):
    """
    :param fastafile: full path to starting FASTA formatted Ensembl peptide file
    :return: None, but create a new identifical FASTA file with the headers updated with location
             information and alternate IDs wherever possible
    """

    # mapping from Ensembl protein ID -> alternate IDs and location
    to_hgnc_id, to_hugo_symbol, to_refseq, to_entrez, to_gene, to_transcript, to_location = ensembl_to_altids()

    # open the new FASTA file to write to
    updated_fastafile = fastafile.replace('.fa', '.withgenelocs.fa')
    out_handle = gzip.open(updated_fastafile, 'w') if updated_fastafile.endswith('.gz') else open(updated_fastafile,
                                                                                                  'w')
    # go through original FASTA file, making updates as needed
    header = ''
    sequence = []

    origfile_handle = gzip.open(fastafile) if fastafile.endswith('gz') else open(fastafile)
    for fasta_line in origfile_handle:
        if fasta_line.startswith('>'):
            if len(sequence) > 0:
                out_handle.write(header + '\n' + ''.join(sequence) + '\n')  # remove all new lines!
            peptide_id = fasta_line[1:-1].split()[0]
            description = {d[:d.find(':')]: d[d.find(':') + 1:] for d in fasta_line[1:-1].split() if ':' in d}
            curr_gene_id = description.get('gene', to_gene.get(peptide_id.split('.')[0], 'UNKNOWN'))

            # scaffolds don't necessarily start with "chromosome"
            chrom_type = [chrom_type for chrom_type, chrom_val in description.items() if BUILD in chrom_val][0]

            sequence = []
            header = ' '.join(['>' + peptide_id,
                               'pep:' + description.get('pep', 'unknown'),
                               chrom_type + ':' + to_location.get(peptide_id.split('.')[0],
                                                               description.get(chrom_type, 'UNKNOWN')),
                               'gene:' + curr_gene_id,
                               'transcript:' + description.get('transcript',
                                                               to_transcript.get(peptide_id.split('.')[0], 'UNKNOWN')),
                               'gene_biotype:' + description.get('gene_biotype', 'UNKNOWN'),
                               'transcript_biotype:' + description.get('transcript_biotype', 'UNKNOWN'),
                               'hgncID:' + to_hgnc_id.get(peptide_id.split('.')[0],
                                                          to_hgnc_id.get(curr_gene_id.split('.')[0], 'UNKNOWN')),
                               'hugoSymbol:' + to_hugo_symbol.get(peptide_id.split('.')[0],
                                                                  to_hugo_symbol.get(curr_gene_id.split('.')[0],
                                                                                     'UNKNOWN'))] +

                              (['refseq:' + to_refseq[peptide_id.split('.')[0]]]
                               if peptide_id.split('.')[0] in to_refseq else []) +
                              (['entrez:' + to_entrez[peptide_id.split('.')[0]]]
                               if peptide_id.split('.')[0] in to_entrez else []))
        else:
            sequence.append(fasta_line.strip())
    origfile_handle.close()

    # Write out the very last sequence (if there is anything to write out)
    if len(sequence) > 0:
        out_handle.write(header + '\n' + ''.join(sequence) + '\n')
        out_handle.close()

    sys.stderr.write('Wrote updated FASTA file to ' + updated_fastafile + '\n')


####################################################################################################

def inflate_toplevel_dnafile(dnafile=data_path + 'ensembl/Homo_sapiens.' + BUILD + '/' +
                                     'Homo_sapiens.' + BUILD + '.dna_sm.toplevel.fa.gz'):
    """
    :param dnafile: full path to a top-level DNA Ensembl file downloaded from, e.g., :
                    ftp://ftp.ensembl.org/pub/release-89/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.toplevel.fa.gz
    :return: None, but inflate the top-level DNA Ensembl file into separate files separated by chromosome
    """

    if not os.path.isfile(dnafile):
        sys.stderr.write('Could not open ' + dnafile + '\n')
        return

    dna_outdir = '/'.join(dnafile.split('/')[:-1])+'/dna_sm/'
    if not os.path.isdir(dna_outdir):
        call(['mkdir', dna_outdir])

    # Needed to keep track of
    chromosome_id = ''
    current_outfile = ''
    chrom_outhandle = None

    dna_inhandle = gzip.open(dnafile) if dnafile.endswith('gz') else open(dnafile)
    for fasta_line in dna_inhandle:
        if fasta_line.startswith('>'):

            # close the current file (finished processing this chromosome)
            if chromosome_id != '':
                chrom_outhandle.write('\n')
                chrom_outhandle.close()
                sys.stderr.write('Wrote to ' + current_outfile + '\n')

            # start a new file for the new chromosome
            chromosome_id = fasta_line.strip().split()[0][1:]
            current_outfile = dna_outdir+'Homo_sapiens.' + BUILD + '.dna_sm.' + chromosome_id + '.fa.gz'
            chrom_outhandle = gzip.open(current_outfile, 'w')
            chrom_outhandle.write(fasta_line)

        else:
            chrom_outhandle.write(fasta_line.strip())  # remove newlines
    dna_inhandle.close()

    # finish the very last file:
    chrom_outhandle.write('\n')
    chrom_outhandle.close()
    sys.stderr.write('Wrote to ' + current_outfile + '\n')


####################################################################################################
# HELPER FUNCTIONS
####################################################################################################

def sequence_ids_by_chromosome(fastafile, chromosome):
    """
    :param fastafile: full path to a FASTA file where location information is found in the header
    :param chromosome: string corresponding to a chromosome (or scaffold/contig) identifier
    :return: set of sequence IDs found on the specified chromosome
    """

    if not os.path.isfile(fastafile):
        sys.stderr.write('Could not open ' + fastafile + '\n')
        return

    sequence_ids = set()  # store protein IDs and transcript IDs
    fasta_inhandle = gzip.open(fastafile) if fastafile.endswith('gz') else open(fastafile)

    for fasta_line in fasta_inhandle:
        if fasta_line.startswith('>'):
            current_chromosome = fasta_line[fasta_line.find('chromosome:'):].split()[0].split(':')[2] \
                if 'chromosome' in fasta_line else None
            if current_chromosome == chromosome:
                sequence_ids.add(fasta_line[1:-1].split()[0])  # add the current protein sequence ID
                if 'transcript:' in fasta_line:
                    sequence_ids.add(fasta_line[fasta_line.find('transcript:') + 11:].split()[0])  # and the cDNA ID
    fasta_inhandle.close()

    return sequence_ids


####################################################################################################

def get_fasta_sequences(fastafile, include_seqs=None):
    """
    :param fastafile: full path to a fasta formatted sequence file
    :param include_seqs: a subset of sequence IDs to return sequences for
    :return: a dictionary of sequenceID->sequence (no newlines)
    """

    if not os.path.isfile(fastafile):
        sys.stderr.write('Could not open ' + fastafile + '\n')
        return

    fasta_inhandle = gzip.open(fastafile) if fastafile.endswith('gz') else open(fastafile)

    seqid_to_sequence = {}  # sequence ID -> full sequence (no newlines)

    sequence_id = ''  # current sequence ID
    sequence = []  # all lines (separated by newline) corresponding to pdb_id

    for fasta_line in fasta_inhandle:
        if fasta_line.startswith('>'):

            # keep track of current stored contents
            if len(sequence) > 0 and (not include_seqs or sequence_id in include_seqs):
                seqid_to_sequence[sequence_id] = ''.join(sequence)

            # reset variables to loop through file:
            sequence_id = fasta_line[1:-1].split()[0]  # remove starting '>' and newline
            sequence = []

        # keep track of the sequence if necessary
        elif not include_seqs or sequence_id in include_seqs:
            sequence.append(fasta_line.strip())
    fasta_inhandle.close()

    # write out remaining contents once we reach the end of the file
    if len(sequence) > 0 and (not include_seqs or sequence_id in include_seqs):
        seqid_to_sequence[sequence_id] = ''.join(sequence)

    return seqid_to_sequence


####################################################################################################

def get_sequence_exons(chrom_seq, location):
    """
    :param chrom_seq: full chromosome sequence (corresponding to the location provided)
    :param location: string containing location information (e.g., "join(123..456,789..1234)" )
    :return: stitch together the exons as specified, returning the DNA sequence, RNA sequence
             and protein sequence. NOTE: Negative indices must be skipped. Exons are already
             in the correct order in the file (and "reversed" if a complementary strand)
    """

    # check if complement is needed
    complement = 'complement' in location
    if complement:
        location = location.replace("complement(", "")

    # exons specified by [start index]..[end index]
    location = location.replace("join(", "").replace(")", "").split(",")

    complete_rna_seq = []
    complete_dna_seq = collections.OrderedDict()

    for exon in location:
        if '..' in exon:
            beg, end = exon.split("..")
            beg = int(NONDECIMAL.sub('', beg))  # exon could start with a "-", dealt with below
            end = int(NONDECIMAL.sub('', end))
        else:  # single nucleotide, not a range
            beg = end = int(NONDECIMAL.sub('', exon))

        # Could be faulty coordinates. Skip if so.
        try:
            exon_subsequence = chrom_seq[beg - 1:end]
        except IndexError:
            exon_subsequence = ''

        if exon.startswith('-'):  # e.g., -3..-3, meaning we need 3 nucleotides here!
            complete_rna_seq.append('N' * beg)
            complete_dna_seq[str(beg) + ':' + str(end)] = 'N' * beg

        elif complement:
            # reverse the exon and replace nucleotides with complementary bases:
            complete_rna_seq.append(''.join([COMPLEMENT[base] for base in exon_subsequence][::-1]))
            complete_dna_seq[str(end) + ':' + str(beg)] = exon_subsequence[::-1]  # reverse indices and exon sequence

        else:
            complete_rna_seq.append(exon_subsequence)
            complete_dna_seq[str(beg) + ':' + str(end)] = exon_subsequence

    complete_rna_seq = ''.join(complete_rna_seq)
    complete_protein_seq = translate_rna(complete_rna_seq.upper())

    return complete_dna_seq, complete_rna_seq, complete_protein_seq


####################################################################################################
# PARSE LOCATION INFORMATION
####################################################################################################

def exons_to_location(exon_list):
    """
    :param exon_list: OrderedDict of exons (as returned by get_sequence_exons)
    :return: string corresponding to those same exons (in the standard CDS format)
    """

    exons = sorted([(int(exon.split(':')[0]), int(exon.split(':')[1])) for exon in exon_list.keys()])
    beg_index, end_index = 0, 1

    # determine if the isoform is on the complementary strand:
    complement = False
    for start, end in exons:
        if start > end:
            complement = True
            break

    if complement:  # reverse the list, so we go LARGEST to SMALLEST index
        exons = exons[::-1]
        beg_index = 1
        end_index = 0

    # now, in either case, add the exons in order
    allexons = ','.join([str(exon[beg_index]) + '..' + str(exon[end_index]) for exon in exons])
    if ',' in allexons:
        allexons = 'join(' + allexons + ')'

    if complement:
        allexons = 'complement(' + allexons + ')'

    # return the string
    return allexons


####################################################################################################

def protein_isoform_location(infile=data_path + 'ensembl/Homo_sapiens.' + BUILD + '/' +
                                    'Homo_sapiens.' + BUILD +'.pep.all.withgenelocs.fa.gz',
                             chrom_id=None):
    """
    :param infile: full path to a FASTA file where the location of each sequence is found in the header
    :param chrom_id: a specific chromosome ID to return results for; if None, return all
    :return: dictionary of chromosome ID -> Ensembl Gene ID -> Ensembl Protein ID ->
             (Ensembl Transcript ID, location, complete header)
    """

    locations = {}  # chromID -> geneID -> protID -> (transID, location, complete header)

    fasta_inhandle = gzip.open(infile) if infile.endswith('gz') else open(infile)
    for fasta_line in fasta_inhandle:
        if fasta_line.startswith('>'):
            ensembl_prot = fasta_line[1:-1].split()[0]  # e.g., ENSP00000263100
            ensembl_gene = fasta_line[fasta_line.find('gene:') + 5:].strip().split()[0]  # e.g., ENSG00000121410
            ensembl_cdna = fasta_line[fasta_line.find('transcript:') + 11:].strip().split()[0]  # e.g., ENST00000263100

            # e.g., GRCh37, 19, complement(join(58864770..58864658,58864658..58864693))
            try:
                (build,
                 chromosome,
                 instructions) = fasta_line[fasta_line.find(BUILD+':'):-1].split()[0].split(':')[:3]
            except ValueError:
                print fasta_line
                sys.exit(1)

            # restrict to a subset of protein IDs if specified
            if chrom_id and chromosome != chrom_id:
                continue

            if chromosome not in locations:
                locations[chromosome] = {}
            if ensembl_gene not in locations[chromosome]:
                locations[chromosome][ensembl_gene] = {}

            locations[chromosome][ensembl_gene][ensembl_prot] = (ensembl_cdna, instructions, fasta_line.strip())
    fasta_inhandle.close()

    # chromID -> geneID -> protID -> (transID, location, complete header)
    return locations


####################################################################################################

def update_sequence_exons(exon_list, rna_alignment):
    """
    :param exon_list: OrderedDict of exons, as returned by get_sequence_exons
    :param rna_alignment: an "alignment" as computed by make_nucleotide_alignment from which exons should be extracted
    :return: updated OrderedDict of exons, where off by 1, 2, or 3 errors in exon start/end positions have been fixed
    """

    # Determine if any exons are reversed
    complement = False
    for exon_range, _ in exon_list.iteritems():
        beg = int(exon_range.split(":")[0])
        end = int(exon_range.split(":")[1])
        if beg > end:
            complement = True
            break

    # if so, get the complement (note that the exon sequence will have already been reversed)
    if complement:
        rna_alignment = ''.join([COMPLEMENT.get(base, base) for base in rna_alignment])

    # Now, build up the new, fixed cds
    fixed_exon_list = collections.OrderedDict()
    i = 0  # exon index
    key = -1
    for exon_range, _ in exon_list.iteritems():
        beg = int(exon_range.split(":")[0])
        end = int(exon_range.split(":")[1])
        if beg > end:
            diff = -1  # how to determine the end of the exon (should we subtract the length, or add it?)
        else:
            diff = 1

        current_exon = ""

        while i < len(rna_alignment) and len(current_exon) <= abs(end - beg):

            # take care of the insertions/deletions (add new exon of the appropriate length!)
            while i < len(rna_alignment) and rna_alignment[i] == '[':

                # first print out what is already sitting in current exon
                if current_exon != "":
                    fixed_exon_list[str(beg) + ':' + str(beg + (len(current_exon) - 1) * diff)] = current_exon
                    beg = beg + len(current_exon) * diff
                    current_exon = ""

                if rna_alignment[i + 1].islower():  # what is inside the []?
                    if str(key) + ':' + str(key) in fixed_exon_list:
                        fixed_exon_list[str(key) + ':' + str(key)] += rna_alignment[i + 1].upper()
                    else:
                        fixed_exon_list[str(key) + ':' + str(key)] = rna_alignment[i + 1].upper()
                else:
                    if current_exon != "":
                        fixed_exon_list[str(beg) + ':' + str(beg + (len(current_exon) - 1) * diff)] = current_exon
                        current_exon = ""
                    beg = beg + (len(current_exon) + 1) * diff
                i = i + 3

            key -= 1

            if i < len(rna_alignment):  # we are still building up this exon!
                current_exon += rna_alignment[i].upper()
                i += 1

        # exited the while loop, is there a new exon to write out!?
        if current_exon != "":
            fixed_exon_list[str(beg) + ':' + str(beg + (len(current_exon) - 1) * diff)] = current_exon

    # we never actually reached the end of the sequence, so make this change...??
    if i < len(rna_alignment):
        fixed_exon_list[str(key) + ':' + str(key)] = rna_alignment[i:].upper().replace('[', '').replace(']', '')

    return fixed_exon_list


####################################################################################################

def translate_rna(rna_seq):
    """
    :param rna_seq: string corresponding to a cDNA sequence
    :return: translation of the cDNA sequence to amino acid sequence using codon_to_aa. Unknown codons or
             accidental stop codons are translated as 'X'
    """

    codon_to_aa = {
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

    prot_seq = []  # amino acids
    for i in range(0, len(rna_seq), 3):

        codon = rna_seq[i:i + 3]
        if len(codon) == 3:
            if codon_to_aa.get(codon, 'X') == "_" and i < len(rna_seq) - 5:  # stop codon somewhere in the middle..
                prot_seq.append('X')
            else:
                prot_seq.append(codon_to_aa.get(codon, 'X'))  # degenerate 3rd positions translate properly

    return ''.join(prot_seq)


####################################################################################################

def transcribe_dna(dna_seq):
    """
    :param dna_seq: OrderedDict of exons (start:end) -> DNA sequence to be stitched together into a single cDNA sequence
    :return: string corresponding to a single concatenated cDNA sequence
    """

    # check if any of the exons are "backwards", as this is an indication of a complement
    complement = False
    for exon_range, exon in dna_seq.iteritems():
        beg = int(exon_range.split(":")[0])
        end = int(exon_range.split(":")[1])
        if beg > end:
            complement = True
            break

    # create a new RNA sequence
    rna_seq = []
    for _, exon in dna_seq.iteritems():
        if complement:  # the exon will already have been reversed! we just need the complementary strand
            rna_seq.append(''.join([COMPLEMENT[i] for i in list(exon)]))
        else:
            rna_seq.append(exon)

    return ''.join(rna_seq)


####################################################################################################

def match_sequences(prot_seq1, prot_seq2, max_miss_rate=MAX_MISS_RATE):
    """
    :param prot_seq1: string corresponding to an amino acid sequence
    :param prot_seq2: string corresponding to an amino acid sequence
    :param max_miss_rate: the maximum miss rate (hamming distance / length) allowed between the two protein
                          sequences to consider them a "good enough" match; anything larger is indicative of
                          something drastically incorrect (frameshift, etc.)
    :return: boolean indicating if the two input protein sequences match well enough
    """

    # only compare strings that are the same length; otherwise they definitely don't match
    if len(prot_seq1) != len(prot_seq2):
        return False

    hamming_distance = 0  # count of mismatches (including gaps) between two sequences
    for i in xrange(len(prot_seq1)):
        if prot_seq1[i] != prot_seq2[i]:
            hamming_distance += 1

    # determine if the proportion of mismatches is within our limit
    return float(hamming_distance) / len(prot_seq1) <= max_miss_rate


####################################################################################################
# ALIGNMENT FUNCTIONS
####################################################################################################

def make_nucleotide_alignment(x, y):
    """
    :param x: nucleotide sequence (possibly containing uppercase and lowercase letters)
    :param y: nucleotide sequence (possibly containing uppercase and lowercase letters)
    :return: an alignment of the two sequences, where deletions are marked by [] and insertions by [x]
    """

    # alignment MATRIX (+1 to avoid IndexError exceptions)
    mat = [[0] * (len(y) + 1) for _ in xrange(len(x) + 1)]  # this keeps track of *scores only*

    # fill in the entire alignment matrix with scores
    for xi in xrange(1, len(x) + 1):
        # it's possible to go "too far" in the x-direction based on the internal y-loop
        diff = 0 if xi == len(x) else 1

        for yi in xrange(1, len(y) + 1):
            mat[xi][yi] = max(mat[xi - 1][yi - 1] + (2 if x[xi - 1].upper() == y[yi - 1].upper() else -1),
                              # allow for a mismatch
                              mat[xi - 1][yi] - 1,  # incorporate a gap
                              mat[xi][yi - 1] - diff,  # incorporate a gap in the other direction
                              0)  # never go below 0

    # indices of where to start backtracking through the matrix:
    xi = 0
    yi = len(y)

    # keep track of the optimal score (with respect to a *GLOBAL* alignment of sequence y)
    best_score = 0
    for i in xrange(len(x) + 1):
        if mat[i][len(y)] > best_score:
            best_score = mat[i][len(y)]
            xi = i

    # backtrack through matrix, filling in the identity of the nucleotide sequence and marking insertions
    best_aligned_seq = []

    while mat[xi][yi] != 0:
        # again, it's possible to go "too far" in the x-direction based on the internal y-loop
        diff = 0 if xi == len(x) else 1

        if mat[xi][yi] == mat[xi - 1][yi] - 1:  # match! leave it capitalized (or as it was)
            best_aligned_seq.insert(0, '[' + x[xi - 1] + ']')
            xi -= 1
        elif mat[xi][yi] == mat[xi][yi - 1] - diff:  # insertion from second sequence! (lowercase)
            best_aligned_seq.insert(0, '[' + y[yi - 1].lower() + ']')
            yi -= 1
        elif mat[xi][yi] == mat[xi - 1][yi - 1] - 1:  # gap in second sequence (lowercase)
            best_aligned_seq.insert(0, x[xi - 1].lower())
            yi -= 1
            xi -= 1
        else:
            best_aligned_seq.insert(0, x[xi - 1])
            yi -= 1
            xi -= 1

    # in case we ended one sequence prematurely, fill in the remaining (non-square matrix)
    while yi != 0 and xi != 0:
        best_aligned_seq.insert(0, x[xi - 1].lower())
        yi -= 1
        xi -= 1
    while yi != 0:
        best_aligned_seq.insert(0, '[' + y[yi - 1].lower() + ']')
        yi -= 1
    while xi != 0:
        best_aligned_seq.insert(0, '[' + x[xi - 1] + ']')
        xi -= 1

    return ''.join(best_aligned_seq), best_score  # return the alignment AND the score


####################################################################################################

def verify_sequence(build, chrom_id, chrom_seq, curr_gene_id, prot_info, all_prot_seqs, all_rna_seqs,
                    exon_dir=data_path + 'ensembl/Homo_sapiens.' + BUILD + '/exons/'):
    """
    :param build: which genome build (e.g., GRCh37, GRCh38)
    :param chrom_id: chromosome ID where the gene is located
    :param chrom_seq: full sequence (no newlines) of the chromosome in nucleotides
    :param curr_gene_id: gene ID that is currently being processed
    :param prot_info: dictionary of protein ID -> (cDNA ID, location, complete header)
    :param all_prot_seqs: dictionary of protein ID -> sequence
    :param all_rna_seqs: dictionary of transcript (cDNA) ID -> sequence
    :param exon_dir: full path to the directory where the output files should be written
    :return: write the protein sequences of the gene to file and return the isoforms that were successfully written
    """

    # keep track of failed protein isoforms:
    progress_handle = gzip.open(PROGRESS_FILE, 'a') if PROGRESS_FILE.endswith('gz') else open(PROGRESS_FILE, 'a')

    for prot_id in sorted(prot_info.keys()):
        # if there is no available protein sequence, move on

        if prot_id in all_prot_seqs:

            # sorted list of exon locs and seqs, full cdna chrom_seq, full prot chrom_seq
            dna_seq, rna_seq, prot_seq = get_sequence_exons(chrom_seq, prot_info[prot_id][1])

            # trim stop codon
            if prot_seq.endswith('_'):
                prot_seq = prot_seq[:-1]

            best_rna_id = prot_info[prot_id][0]  # best cDNA ID

            fix = False

            # if the translated protein doesn't match what it should have...
            if not match_sequences(prot_seq.upper(), all_prot_seqs[prot_id].upper(), MAX_MISS_RATE):

                best_alignment_score = 0  # best alignment score (2 for a match, -1 for anything else)
                best_rna_alignment = ""  # error checking (print "alignment" to log file) and updating
                expected_rna_length = 3 * len(all_prot_seqs[prot_id])  # how long *should* the transcript be?

                # go through all other transcript IDs in case we were considering the wrong one
                for curr_rna_id, curr_rna_seq in all_rna_seqs.items():

                    # and also allow for a shift (if possible)
                    for s in xrange(len(curr_rna_seq) - expected_rna_length + 1):

                        # get the protein sequence from translating this (possibly shifted) cDNA sequence
                        # NOTE: by definition, this will be the same length as the expected protein
                        curr_prot_seq = translate_rna(curr_rna_seq[s:s + expected_rna_length].upper())

                        # in the event that the shifted transcript now matches the protein as it should:
                        if match_sequences(curr_prot_seq.upper(), all_prot_seqs[prot_id].upper(), MAX_MISS_RATE):

                            # align the transcript sequence here to what we should have gotten from the DNA; usually
                            # these two sequences are off by up to 2bp (frameshift), which is why we +3 at the end:
                            (curr_rna_alignment,
                             current_alignment_score) = make_nucleotide_alignment(
                                rna_seq, curr_rna_seq[s:s + expected_rna_length + 3])
                            if current_alignment_score > best_alignment_score:
                                best_rna_alignment = curr_rna_alignment
                                best_alignment_score = current_alignment_score
                                best_rna_id = curr_rna_id

                # Was our "best fix" within the maximum miss rate?
                # NOTE: the alignment scores give a value of '2' for every match (thus the *2)
                if float(best_alignment_score) >= (1 - MAX_MISS_RATE) * 2 * len(rna_seq):
                    dna_seq = update_sequence_exons(dna_seq, best_rna_alignment)  # new ordered list of fixed exons
                    rna_seq = transcribe_dna(dna_seq)
                    prot_seq = translate_rna(rna_seq.upper())

                    if prot_seq.endswith('_'):  # trim stop codon
                        prot_seq = prot_seq[:-1]

                    fix = True

                if not fix:  # this particular gene, transcript, protein failed: print out the information:
                    progress_handle.write(curr_gene_id + ', ' + best_rna_id + ', ' + prot_id + '\t' +
                                          'ALN_SCORE=' + str(best_alignment_score) + '/' +
                                          str(2 * len(rna_seq)) + '\t' +
                                          'PASS_STATUS=' + str(fix) + '\n' +
                                          'ALN=' + best_rna_alignment + '\n')
            else:  # proteins already matched!
                fix = True

            # if it already matched OR we fixed it
            if fix:
                # write out, in fasta format, the corresponding protein and cDNA files:
                for sequence, sequence_type in [(prot_seq, 'prot'), (rna_seq, 'cdna')]:
                    sequence_handle = open(exon_dir + chrom_id + '/' + curr_gene_id + '/' +
                                           prot_id + '.' + sequence_type + '.fa', 'w')
                    sequence_handle.write('>' + prot_id + ' prot:' + prot_id + ' ')

                    # e.g., gene:ENSGxx refseq:X hugoID:X
                    descriptors = prot_info[prot_id][2].split()[1:] + ['length:' + str(len(sequence))]
                    sequence_handle.write(
                        ' '.join(['chromosome:' + build + ':' + chrom_id + ':' + exons_to_location(dna_seq)
                                  if d.startswith('chromosome') else d for d in descriptors]) +
                        '\n' + sequence + '\n\n')
                    sequence_handle.close()

                # write out the exons, in order, for easier translation of mutations (in the end):
                exon_handle = open(exon_dir + chrom_id + '/' + curr_gene_id + '/' + prot_id + '.exons.txt', 'w')
                exon_handle.write('>' + prot_id + ' prot:' + prot_id + ' ')
                descriptors = prot_info[prot_id][2].split()[1:] + ['length:' + str(len(rna_seq))]
                exon_handle.write(' '.join(['chromosome:' + build + ':' + chrom_id + ':' + exons_to_location(dna_seq)
                                            if d.startswith('chromosome') else d for d in descriptors]) + '\n')
                for exon_range, exon in dna_seq.iteritems():
                    exon_handle.write(exon_range + '\t' + exon + '\n')
                exon_handle.close()

    progress_handle.close()


####################################################################################################
# POST-PROCESSING: slim down peptide file for subsequent analyses
####################################################################################################

def create_nonredundant_sequences(infile):
    """
    :param infile: full path to a fasta file containing potentially redundant sequences
    :return: None, but write out to a new file a nonredundant version of the original file
    """

    # process the original sequences, storing all sequences and their corresponding identifiers
    all_sequences = {}  # original sequence -> set(identifiers)

    orig_fasta = gzip.open(infile) if infile.endswith('gz') else open(infile)

    sequence_id, current_seq = '', []
    for fasta_line in orig_fasta:

        if fasta_line.startswith('>'):
            if len(current_seq) > 0:
                sequence = ''.join(current_seq)
                if sequence not in all_sequences:
                    all_sequences[sequence] = set()
                all_sequences[sequence].add(sequence_id)

            # and now, reset:
            sequence_id = fasta_line[1:-1].split()[0]  # original sequence identifier
            current_seq = []

        else:
            current_seq.append(fasta_line.strip().upper().replace('*', 'X').replace('_', ''))
    orig_fasta.close()

    # very last sequence?
    if len(current_seq) > 0:
        sequence = ''.join(current_seq)
        if sequence not in all_sequences:
            all_sequences[sequence] = set()
        all_sequences[sequence].add(sequence_id)

    # write out a new file containing only nonredundant sequences
    unique_sequences = len(str(len(all_sequences.keys())))  # length of longest new sequence ID

    nr_fasta = infile.replace('.fa', '.nonredundant.fa')  # filename for the nonredundant version of the FASTA file
    out_handle = gzip.open(nr_fasta, 'w') if nr_fasta.endswith('gz') else open(nr_fasta, 'w')
    for new_seq_id, full_sequence in enumerate(sorted(all_sequences.keys())):
        if len(full_sequence) > 0:  # only bother with sequences with at least one amino acid...
            out_handle.write('>seq' + str(new_seq_id).zfill(unique_sequences) + ' ' +
                             ','.join(sorted(list(all_sequences[full_sequence]))) + '\n' +
                             full_sequence + '\n\n')
    out_handle.close()

    sys.stderr.write('Nonredundant version of ' + infile + ' found in ' + nr_fasta + '\n')


####################################################################################################
# MAIN
####################################################################################################

if __name__ == "__main__":

    # Parse the command-line arguments
    parser = argparse.ArgumentParser(description='Verify all Ensembl protein sequences using the peptide, cDNA ' +
                                                 'and DNA sequences.')

    parser.add_argument('--chromosome', type=str,
                        help='Chromosome or scaffold/contig to verify protein sequences on.',
                        default='Y')
    parser.add_argument('--ensembl_dir', type=str,
                        help='Directory to where all sequencing information from Ensembl is and will be stored',
                        default=data_path + 'ensembl/Homo_sapiens.' + BUILD + '/')
    parser.add_argument('--fix_fasta', dest='fix_fasta', action='store_true', default=False,
                        help='Add alternate IDs and location information to the original Ensembl protein file.')
    parser.add_argument('--inflate_toplevel', dest='inflate_toplevel', action='store_true', default=False,
                        help='"Inflate" the originally downloaded soft-masked DNA file into separate DNA files '
                             'for each chromosome or scaffold.')
    parser.add_argument('--verify_exons', dest='verify_exons', action='store_true', default=False,
                        help='Create new FASTA files of protein sequences, cDNA sequences, and exons where only '
                             'those protein isoforms whose DNA sequences matched the corresponding cDNA sequences '
                             'and whose cDNA sequences translated to the corresponding protein sequences are '
                             'retained.')
    parser.add_argument('--create_final_fasta', dest='create_final_fasta', action='store_true', default=False,
                        help='Create a new FASTA-formatted protein sequence file containing the subset of protein'
                             'isoforms that we were able to successfully verify.')
    args = parser.parse_args()

    # ----------------------------------------------------------------------------------------------------
    if args.fix_fasta:
        """
        Create a new FASTA file where alternate IDs and location information has been added to the header and sequences
        are found on single lines.
    
        NOTE: We need the following files obtained from Ensembl's BioMart (see function description for more details):
        1) hgncfile=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.toHGNC.tsv', 
        2) refseqfile=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.toRefSeq.tsv',
        3) entrezfile=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.toEntrez.tsv',
        4) locfile=data_path+'ensembl/Homo_sapiens.GRCh38/Homo_sapiens.GRCh38.genelocs.tsv'
        """

        update_fasta_genelocs_altids(args.ensembl_dir + 'Homo_sapiens.' + BUILD + '.pep.all.fa.gz')

    # ----------------------------------------------------------------------------------------------------
    elif args.inflate_toplevel:
        """
        All the chromosome files can be downloaded from Ensembl. The separate files only include chromosomes 
        1-22, X, Y, MT though. In order to get the same fasta files for non-chromosomal DNA, the "top level" file must 
        be downloaded and processed instead.
        """

        inflate_toplevel_dnafile(args.ensembl_dir + 'Homo_sapiens.' + BUILD + '.dna_sm.toplevel.fa.gz')

    # ----------------------------------------------------------------------------------------------------
    elif args.verify_exons:
        """
        With a FASTA file of protein sequences *with gene locations* and alternate IDs, verify all sequences. The 
        soft-masked DNA sequences and cDNA sequences must have been downloaded from Ensembl.
        """

        protein_file = args.ensembl_dir + 'Homo_sapiens.' + BUILD + '.pep.all.withgenelocs.fa.gz'

        # make sure that there are genes to process on the specified chromosome:
        mapping = protein_isoform_location(protein_file, args.chromosome)

        # Nothing to process on this chromosome, exit
        if args.chromosome not in mapping:
            sys.stderr.write('No genes found on chromosome ' + args.chromosome + '\n')
            sys.exit(0)

        # Otherwise, set up the output directory
        for exon_dir in [args.ensembl_dir + 'exons/', args.ensembl_dir + 'exons/' + args.chromosome]:
            if not os.path.isdir(exon_dir):
                call(['mkdir', exon_dir])

        # reset the progress file and make sure that it is blank:
        PROGRESS_FILE = PROGRESS_FILE.replace('.log',
                                              '-' + args.chromosome + '.log').replace('/exons/',
                                                                                      '/exons/'+args.chromosome+'/')
        log_handle = gzip.open(PROGRESS_FILE, 'w') if PROGRESS_FILE.endswith('gz') else open(PROGRESS_FILE, 'w')
        log_handle.close()

        # store the protein sequences and cDNA sequences *for our chromosome of interest*:
        include_sequences = sequence_ids_by_chromosome(protein_file, args.chromosome)

        prot_seqs = get_fasta_sequences(protein_file, include_sequences)
        cdna_seqs = get_fasta_sequences(args.ensembl_dir + 'Homo_sapiens.' + BUILD + '.cds.all.fa.gz',
                                        include_sequences)

        # store the entire chromosome (DNA sequence). There should be only one sequence per file, but we double check
        # (just in case)
        dna_file_path = args.ensembl_dir + 'dna_sm/Homo_sapiens.' + BUILD + '.dna_sm.' + args.chromosome + '.fa.gz'

        chromosome_sequence = get_fasta_sequences(dna_file_path, [args.chromosome])[args.chromosome]

        # for each gene on this particular chromosome
        for gene_id in sorted(mapping[args.chromosome]):

            if not os.path.isdir(args.ensembl_dir + 'exons/' + args.chromosome + '/' + gene_id):
                call(['mkdir', args.ensembl_dir + 'exons/' + args.chromosome + '/' + gene_id])

            curr_prot_seqs = {pid: prot_seqs[pid].replace('*', 'X') for pid in mapping[args.chromosome][gene_id].keys()}
            curr_cdna_seqs = {tid: cdna_seqs[tid] for tid in [a[0] for a in mapping[args.chromosome][gene_id].values()]}

            sys.stderr.write(args.chromosome + ', ' + gene_id + '...')
            verify_sequence(BUILD, args.chromosome, chromosome_sequence, gene_id, mapping[args.chromosome][gene_id],
                            curr_prot_seqs, curr_cdna_seqs, args.ensembl_dir + 'exons/')
            sys.stderr.write('Done.\n')

        sys.stderr.write('Logfile in ' + PROGRESS_FILE + '\n')

    # ----------------------------------------------------------------------------------------------------
    elif args.create_final_fasta:
        """
        Create a "verified" FASTA protein file as well as a nonredundant version of this file (to speed up subsequent 
        analyses)
        """

        original_fasta = args.ensembl_dir + 'Homo_sapiens.' + BUILD + '.pep.all.withgenelocs.fa.gz'
        new_prot_fasta = original_fasta.replace('.fa', '.verified.fa')
        new_cdna_fasta = original_fasta.replace('.fa', '.verified.fa').replace('.pep.', '.cdna.')

        sys.stderr.write('Attempting to create '+new_prot_fasta+'\n')
        sys.stderr.write('Attempting to create '+new_cdna_fasta+'\n')

        chromosomes = [str(chrom_id) for chrom_id in range(1, 23)] + ['X', 'Y', 'MT']
        alt_chroms = sorted([dir_name for dir_name in os.listdir(args.ensembl_dir+'exons/')
                             if os.path.isdir(args.ensembl_dir+'exons/'+dir_name) and
                             dir_name not in chromosomes])

        prot_fa_handle = gzip.open(new_prot_fasta, 'w') if new_prot_fasta.endswith('gz') else open(new_prot_fasta, 'w')
        cdna_fa_handle = gzip.open(new_cdna_fasta, 'w') if new_cdna_fasta.endswith('gz') else open(new_cdna_fasta, 'w')
        for chrom_id in chromosomes+alt_chroms:
            for gene_id in sorted([dir_name for dir_name in os.listdir(args.ensembl_dir+'exons/'+chrom_id)
                                   if os.path.isdir(args.ensembl_dir+'exons/'+chrom_id+'/'+dir_name)]):

                # keep track of protein sequences:
                for prot_file in sorted([file_name for file_name in
                                         os.listdir(args.ensembl_dir+'exons/'+chrom_id+'/'+gene_id)
                                         if file_name.endswith('.prot.fa')]):
                    with open(args.ensembl_dir+'exons/'+chrom_id+'/'+gene_id+'/'+prot_file) as prot_handle:
                        for prot_line in prot_handle:
                            prot_fa_handle.write(prot_line)

                # and also newly verified cdna sequences:
                for cdna_file in sorted([file_name for file_name in
                                         os.listdir(args.ensembl_dir+'exons/'+chrom_id+'/'+gene_id)
                                         if file_name.endswith('.cdna.fa')]):
                    with open(args.ensembl_dir+'exons/'+chrom_id+'/'+gene_id+'/'+cdna_file) as cdna_handle:
                        for cdna_line in cdna_handle:
                            cdna_fa_handle.write(cdna_line)
        prot_fa_handle.close()
        cdna_fa_handle.close()

        sys.stderr.write('Verified protein sequences in ' + new_prot_fasta + '\n')
        sys.stderr.write('Verified cDNA sequences in ' + new_cdna_fasta + '\n')

        # create_nonredundant_sequences(new_prot_fasta)
