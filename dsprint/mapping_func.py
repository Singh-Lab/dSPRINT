import pandas as pd
import numpy as np
import pickle
import math


def correct_exons_frameshift(exon_df, targetid, frameshift_file):
    """
    Correct the exons end/start positions if there's a frameshift.
    Frameshift are extracted from Shilpa's exons sequences.
    This function corrects the indices of the exons table accordingly.
    """

    # Get the frameshifts index and length of the exons
    with open(frameshift_file, 'rb') as f:
        exons_frameshifts = pickle.load(f)

    for frameshift in exons_frameshifts[targetid + ".exons.txt"]:
        idx, length, _ = frameshift

        # Find the exon we need to add bps to
        first_bp_count = 1
        for index, exon in exon_df.iterrows():
            ex_start = int(exon[0])
            ex_end = int(exon[1])
            exon_len = (ex_end - ex_start + 1)

            # Fixing start pos of the exon
            if idx <= first_bp_count:
                exon_df.loc[index].start_pos = ex_start - length
                break
            # Fixing end pos of the exon
            elif idx <= (first_bp_count + exon_len):
                exon_df.loc[index].end_pos = ex_end + length
                break
            first_bp_count += exon_len


# -------------------------------------------------------------------------------------------#

def create_exon_pos_table(chrom_raw, targetid, frameshift_file):
    """
    A function that get chromosome raw data from the hmmer results and return a data-frame of the exons.
    """
    exons_raw = chrom_raw

    # Removing the complement bracates if exist
    if exons_raw.find("complement(") >= 0:
        exons_raw = exons_raw[exons_raw.find("complement(") + 11:-1]

    # Removing the join bracates if exist
    if exons_raw.find("join(") >= 0:
        exons_raw = exons_raw[exons_raw.find("join(") + 5:-1]

    # In case there's only one exon, take everything after the second ":"
    else:
        exons_raw = exons_raw[exons_raw.find(":", chrom_raw.find(":") + 1) + 1:]

    exons_list = exons_raw.split(",")
    exon_pos = []
    frameshift_flag = False
    for ex in exons_list:
        # flag cases where Shilpa added "-" to a position number to signify frameshift in the sequences
        if ex[0] == "-":
            frameshift_flag = True
            continue

        # Adding the real exons to exons_pos list
        exon_pos.append([int(p) for p in ex.split("..")])

    # Creating a table for the start and end of exons
    exon_df = pd.DataFrame(exon_pos, columns=('start_pos', 'end_pos'))

    # Correct frameshift if frameshift exist
    if frameshift_flag:
        correct_exons_frameshift(exon_df, targetid, frameshift_file)

    exon_df['length'] = exon_df['end_pos'] - exon_df['start_pos'] + 1
    exon_df['first_bp_count'] = exon_df.length.cumsum() - exon_df.length + 1
    return exon_df


# -------------------------------------------------------------------------------------------#

def find_protein_pos(chrom_pos, exon_df, chrom_raw):
    """
    A function that get chromosome position and data-frame of exons,
    and return the protein position or -1 if it's not within any exon.
    Not currently beeing used in the pipeline.
    """
    # Go over all the exons to search for the chrom position there
    for index, exon in exon_df.iterrows():
        start_pos = int(exon[0])
        end_pos = int(exon[1])
        first_bp_count = int(exon[3])

        # If the chrom position is inside this exon
        if start_pos <= chrom_pos <= end_pos:

            # Calculate position for reverse complement strand:
            # the protein is translated from the end position towards the start position of the exon
            if chrom_raw.find("complement") >= 0:
                len_from_exon_start = end_pos - chrom_pos
            # Calculate position for forward strand
            else:
                len_from_exon_start = chrom_pos - start_pos

            # Calculate the position on the mRNA transcript
            transcript_pos = len_from_exon_start + first_bp_count

            # Calculate the position on the protein sequence
            protein_pos = int(math.ceil(float(transcript_pos) / 3))

            return protein_pos

    # If the position wasn't in the regions of any exon
    return -1


def find_chrom_bps(protein_pos, exon_table, chrom_raw_data):
    """
    A function that gets protein position and data-frame of exons,
    and returns the chromosome positions of the corresponding codon, as well as the nucleobases
    at those positions
    """
    seq = ''

    # calculate the mRNA transcript index of this protein position (the 1st bp in the triplet)
    transcript_pos = (protein_pos * 3) - 2

    # Iterating over all the gene exons
    for index, exon in exon_table.iterrows():
        first_bp_count = int(exon["first_bp_count"])
        exon_length = int(exon["length"])
        last_bp_count = first_bp_count + exon_length - 1

        # Checking if the transcript position is within this exon
        if first_bp_count <= transcript_pos <= last_bp_count:

            start_pos = int(exon["start_pos"])
            end_pos = int(exon["end_pos"])

            len_from_exon_start = transcript_pos - first_bp_count

            # Calculate bps position for reverse complement strand:
            # the protein is translated from the end position towards the start position of the exon
            if chrom_raw_data.find("complement") >= 0:
                chrom_pos_1st = end_pos - len_from_exon_start
                seq += exon.seq[chrom_pos_1st - start_pos]

                chrom_pos_2nd = chrom_pos_1st - 1
                # If the exons end here: move to the next exon
                if chrom_pos_2nd < start_pos:
                    index += 1
                    chrom_pos_2nd = int(exon_table["end_pos"][index])
                    start_pos = int(exon_table["start_pos"][index])

                seq += exon_table["seq"][index][chrom_pos_2nd - start_pos]

                # If the exons ends here: move to the next exon
                chrom_pos_3rd = chrom_pos_2nd - 1
                if chrom_pos_3rd < start_pos:
                    index += 1
                    chrom_pos_3rd = int(exon_table["end_pos"][index])
                    start_pos = int(exon_table["start_pos"][index])

                seq += exon_table["seq"][index][chrom_pos_3rd - start_pos]

            # Calculate position for forward strand
            else:
                chrom_pos_1st = start_pos + len_from_exon_start
                seq += exon.seq[chrom_pos_1st - start_pos]

                chrom_pos_2nd = chrom_pos_1st + 1
                # If the exons end here: move to the next exon
                if chrom_pos_2nd > end_pos:
                    index += 1
                    start_pos = int(exon_table["start_pos"][index])
                    end_pos = int(exon_table["end_pos"][index])
                    chrom_pos_2nd = start_pos

                seq += exon_table["seq"][index][chrom_pos_2nd - start_pos]

                # If the exons end here: move to the next exon
                chrom_pos_3rd = chrom_pos_2nd + 1
                if chrom_pos_3rd > end_pos:
                    index += 1
                    start_pos = int(exon_table["start_pos"][index])
                    chrom_pos_3rd = start_pos

                seq += exon_table["seq"][index][chrom_pos_3rd - start_pos]

            return (chrom_pos_1st, chrom_pos_2nd, chrom_pos_3rd), seq


def is_number(s):
    """
    Boolean function - determine if a given text can be converted to a number
    """
    try:
        float(s)
        return True
    except ValueError:
        return False
