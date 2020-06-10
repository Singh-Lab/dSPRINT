import sys

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from dsprint.core import retrieve_exon_seq


# A dictionary that maps codons to amino acids
codon_table = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}


def create_alt_codon(exac_ref_bp, curr_alt_bp, ref_codon, alt_codon_pos, is_complementary=False):
    """
    A function that create the new codon for the alteration
    """

    # For error logging
    functionNameAsString = sys._getframe().f_code.co_name

    # Complement strand - transversing the bp to base-complement
    if is_complementary:
        new_bp = ""
        for c in curr_alt_bp:
            if (c.upper() == 'A'):
                new_bp = new_bp + 'T'
            elif (c.upper() == 'T'):
                new_bp = new_bp + 'A'
            elif (c.upper() == 'G'):
                new_bp = new_bp + 'C'
            else:
                new_bp = new_bp + 'G'
            new_bp = new_bp[::-1]  # TODO: is the reverse needed?

        exac_ref_bp_adj = ""
        for c in exac_ref_bp:
            if (c.upper() == 'A'):
                exac_ref_bp_adj = exac_ref_bp_adj + 'T'
            elif (c.upper() == 'T'):
                exac_ref_bp_adj = exac_ref_bp_adj + 'A'
            elif (c.upper() == 'G'):
                exac_ref_bp_adj = exac_ref_bp_adj + 'C'
            else:
                exac_ref_bp_adj = exac_ref_bp_adj + 'G'
            exac_ref_bp_adj = exac_ref_bp_adj[::-1]  # TODO: is the reverse needed?

    # Regular strand
    else:
        new_bp = curr_alt_bp
        exac_ref_bp_adj = exac_ref_bp

    # Validation: making sure the ref bp from ExAC is inside the ref codon sequence retrieved from hg19 or the other way around (at least one contain the other)
    if (ref_codon.find(exac_ref_bp_adj) == -1 and exac_ref_bp_adj.find(ref_codon) == -1):
        print
        functionNameAsString + " Error: ExAC ref sequence " + exac_ref_bp_adj + " isn't found in hg19 retrieved codon sequence " + ref_codon

    new_alt_codon = ref_codon[:alt_codon_pos] + new_bp + ref_codon[alt_codon_pos + len(exac_ref_bp_adj):]
    if (len(new_alt_codon) != 3):
        print
        functionNameAsString + " Error: new alt codon length isn't 3: " + new_alt_codon

    return new_alt_codon


# -------------------------------------------------------------------------------------------#

def retrieve_codon_seq(chrom_pos_list, chrom, hg19_file, is_complementary=False):
    """
    Retrieve the codon base-pairs from the ref sequence
    """
    assert(len(chrom_pos_list) == 3)

    # reduce individual calls to retrieve_exon_seq if possible, for speedup
    is_contiguous = chrom_pos_list[0] == chrom_pos_list[1] - 1 == chrom_pos_list[2] - 2
    if is_contiguous:
        seq = retrieve_exon_seq(chrom_pos_list[0], chrom_pos_list[2], chrom, hg19_file)
    else:
        seq = ''.join([retrieve_exon_seq(pos, pos, chrom, hg19_file) for pos in chrom_pos_list])

    if is_complementary:
        return str(Seq(seq, generic_dna).complement())
    else:
        return seq


# -------------------------------------------------------------------------------------------#

def exac_validation_checks(chrom_alter, protein_pos, aa, alt_codon_pos, chrom_pos, bp_ref):
    """
    Validation checks on the data received from ExAC
    :param chrom_alter:
    :param protein_pos:
    :param aa:
    :param alt_codon_pos:
    :param chrom_pos:
    :param bp_ref:
    :return:
    """
    error_flag = False

    # Validation: the ExAC chromosome position is within a protein
    exac_prot_data = True
    if chrom_alter["PROTEIN_POSITION"] == "":
        print(" Error: ExAC chromosome position " + str(chrom_pos) + " doesn't correspond to a protein")
        # We assume it's an error in ExAC and logging alteration anyway.
        exac_prot_data = False

    else:
        # Validation: the ExAC protein position match the HMMER protein position
        exac_prot_pos = chrom_alter["PROTEIN_POSITION"]
        # in case there's more than one position listed
        if exac_prot_pos.find("-") != -1:
            # Trying to convert to int, but leaving as string if its a question mark
            try:
                first_exac_prot_pos = int(exac_prot_pos[:exac_prot_pos.find("-")])
            except:
                first_exac_prot_pos = exac_prot_pos[:exac_prot_pos.find("-")]
            try:
                last_exac_prot_pos = int(exac_prot_pos[exac_prot_pos.find("-") + 1:])
            except:
                last_exac_prot_pos = exac_prot_pos[exac_prot_pos.find("-") + 1:]
        else:
            first_exac_prot_pos = int(float(exac_prot_pos))  # float if for when number is listed in ExAC with dd.0
            last_exac_prot_pos = first_exac_prot_pos
        # Checking of the protein position isn't within the range described by ExAC
        if not (first_exac_prot_pos <= protein_pos <= last_exac_prot_pos):
            print(str(chrom_pos) + " Error: ExAC protein position " + str(
                first_exac_prot_pos) + " doesn't match HMMER protein position " + str(protein_pos))
            error_flag = True

        # Validation: the ExAC aa match the HMMER aa
        exac_aa = chrom_alter["AMINO_ACIDS"]
        if exac_aa.find("/") != -1:
            exac_ref_aa = exac_aa[:exac_aa.find("/")]
        else:
            exac_ref_aa = exac_aa
        exac_alt_aa = exac_aa[exac_aa.find("/") + 1:]
        if exac_ref_aa != aa:
            if exac_ref_aa == "":
                print(str(chrom_pos) + " ExAC amino acid is blank")
            else:
                print(" Error: ExAC amino acid identity " + exac_ref_aa + " doesn't match HMMER amino-acid " + aa)
                error_flag = True

        # Extracting aa codon data if exist
        exac_codons = chrom_alter["CODONS"]
        exac_ref_codon = exac_codons[:exac_codons.find("/")]
        exac_alt_codon = exac_codons[exac_codons.find("/") + 1:]
        # Validation: the ExAC codon match the returned codon sequence from hg19, or at least one contain the other
        if ((exac_ref_codon.upper().find(bp_ref.upper()) == -1) and (
                bp_ref.upper().find(exac_ref_codon.upper()) == -1)):
            print(str(chrom_pos) + " Error: ExAC bp codon " + exac_ref_codon.upper() + " doesn't match hg19 codon sequence retrieved " + bp_ref.upper())

        return exac_prot_data, exac_alt_aa, exac_alt_codon, error_flag