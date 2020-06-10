import sys
import pandas as pd
from enum import Enum


# Constants to represent indels types
class indel_type(Enum):
    FRAME_SHIFT_INDEL = 1
    IN_FRAME_INDEL = 2
    NO_INDEL = 0


# -------------------------------------------------------------------------------------------#

# A boolean function that check according to ref and alt if it's frameshift indel or not
def is_indel(ref, alt, chrom_alter):
    # For error logging
    functionNameAsString = sys._getframe().f_code.co_name

    if (len(ref) != len(alt)):
        if (len(ref) > len(alt)):
            indel_len = len(ref) - len(alt)
        else:
            indel_len = len(ref) - len(alt)
        if ((indel_len % 3) > 0):

            # Validation: ExAC identify this as a frameshift too
            conseq = chrom_alter['CONSEQUENCE']
            # if ("frameshift" not in conseq):
            #   print functionNameAsString+" Error: ExAC doesn't recognize a frameshift indel - "+str(chrom_alter["pos"])

            # frameshift indel
            return indel_type.FRAME_SHIFT_INDEL

        # In frame indel
        else:
            return indel_type.IN_FRAME_INDEL

    # Not an indel
    else:
        return indel_type.NO_INDEL


# -------------------------------------------------------------------------------------------#

def diff(a, b):
    """
    A function that returns the indices of the strings differences
    """
    return [i for i in range(len(a)) if a[i] != b[i]]


def table_editing(chrom_gene_table):
    """
    A function that adds more lines for all the effected chromosomal positions of each indel
    """
    i_cnt = 0
    d_cnt = 0
    indels_table = pd.DataFrame(columns=chrom_gene_table.columns)
    indels_table_i = 0
    comments_col = []

    for index, line in chrom_gene_table.iterrows():
        ref = line["REF"]
        alt = line["ALT"]
        pos = line["POS"]
        strand = line["STRAND"]

        # Handling deletion
        if len(ref) > len(alt):
            d_cnt += 1
            comments_col.append("ignore for this position: d-" + str(d_cnt))
            # Adding to indels table only the inframe ones
            if is_indel(ref, alt, line) == indel_type.IN_FRAME_INDEL:
                init_pos = pos + len(alt)
                deletion_len = len(ref) - len(alt)
                for j in range(deletion_len):
                    new_line = line.copy(deep=True)
                    new_line["POS"] = init_pos + j
                    new_line["REF"] = ref[j + 1]
                    new_line["ALT"] = "-"
                    new_line["COMMENTS"] = "d-" + str(d_cnt)
                    indels_table.loc[indels_table_i] = new_line
                    indels_table_i += 1

        # Handling insertion
        elif len(alt) > len(ref):
            i_cnt += 1
            comments_col.append("ignore for this position: i-" + str(i_cnt))

            # Adding to indels table only the inframe ones
            if is_indel(ref, alt, line) == indel_type.IN_FRAME_INDEL:
                init_pos = pos + len(ref)
                # insertion_len = len(alt) - len(ref)
                new_line = line.copy(deep=True)
                if strand == 1:
                    new_line["POS"] = init_pos
                else:
                    new_line["POS"] = init_pos - 1
                new_line["REF"] = "-"
                new_line["ALT"] = alt[len(ref):]
                new_line["COMMENTS"] = "i-" + str(i_cnt)
                indels_table.loc[indels_table_i] = new_line
                indels_table_i += 1

        # Handling mismatch written with redundant bps
        elif len(ref) > 1:
            diff_idx = diff(ref, alt)
            # A case when only the first bp is the alteration
            if diff_idx == [0]:
                # Fix ref and alt fields
                chrom_gene_table.at[index, "REF"] = ref[0]
                chrom_gene_table.at[index, "ALT"] = alt[0]

                strand = line["STRAND"]

                # Fix amino_acids field
                aa = line["AMINO_ACIDS"]
                # If aa field is not empty
                if aa != "":
                    if strand == 1:
                        if aa.find("/") != -1:
                            new_aa = aa[0] + aa[aa.find("/"):aa.find("/") + 2]
                        else:
                            new_aa = aa[0]
                    else:
                        if aa.find("/") != -1:
                            new_aa = aa[aa.find("/") - 1:aa.find("/") + 1] + aa[-1]
                        else:
                            new_aa = aa[-1]
                    chrom_gene_table.at[index, "AMINO_ACIDS"] = new_aa

                # Fix prot_pos field
                prot_pos = line["PROTEIN_POSITION"]
                if prot_pos != "":
                    if prot_pos.find("-") != -1:
                        if strand == 1:
                            chrom_gene_table.at[index, "PROTEIN_POSITION"] = prot_pos[:prot_pos.find("-")]
                        else:
                            chrom_gene_table.at[index, "PROTEIN_POSITION"] = prot_pos[prot_pos.find("-") + 1:]

                comments_col.append("removed redundant bps")
            else:
                comments_col.append("")

        # No Indel
        else:
            comments_col.append("")

    chrom_gene_table["COMMENTS"] = comments_col
    # print "Number of insertions = "+str(i_cnt)
    # print "number of deletion = "+str(d_cnt)

    return indels_table
