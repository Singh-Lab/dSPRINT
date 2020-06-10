import pandas as pd
from collections import defaultdict

try:
    snakemake
except NameError:
    import sys

    if len(sys.argv) != 3:
        print('Usage: <script> <input_hmmer_csv> <output_hmmer_csv>')
        sys.exit(0)

    INPUT_FILE, OUTPUT_FILE = sys.argv[1:]
else:
    INPUT_FILE = snakemake.input[0]
    OUTPUT_FILE = snakemake.output[0]

if __name__ == '__main__':

    df = pd.read_csv(INPUT_FILE, sep='\t', skiprows=[0, 1], header=0)

    df['pfam_id'], df['domain_name'] = zip(*df['HMM_Name'].apply(lambda x: x.split('_', 1)))
    del df['HMM_Name']
    # Get the columns to the original order
    cols = df.columns.tolist()
    cols = cols[0:1] + cols[-2:] + cols[1:-2]
    df = df[cols]

    # Parsing the description field
    description_headers = ["prot", "pep", "chromosome", "gene", "transcript", "gene_biotype", "transcript_biotype",
                           "hgncID", "hugoSymbol", "refseq", "entrez", "length", "HMMStart", "HMMEnd"]
    descriptions_parsed = defaultdict(list)
    desc_list = df["Description"].tolist()

    # Iterating over the Description field and parsing the data
    for i in range(len(desc_list)):
        curr_desc = desc_list[i].split(' ')
        desc_idx = 0
        # Iterating over the headers list and looking for the headers in the parsed data
        for headers_idx in range(len(description_headers)):
            val = curr_desc[desc_idx]
            if val.find(description_headers[headers_idx]) == 0:
                if description_headers[headers_idx] in ("HMMStart", "HMMEnd"):
                    parsed = val[val.find("=") + 1:-1]
                else:
                    parsed = val[val.find(":") + 1:]
                descriptions_parsed[description_headers[headers_idx]].append(parsed)
                desc_idx += 1
            else:
                # If the headers isn't found, appending empty string instead
                descriptions_parsed[description_headers[headers_idx]].append("")

    # Adding the parsed description to the data frame
    del df["Description"]
    for col_name in description_headers:
        df[col_name] = descriptions_parsed[col_name]

    # Seperate chromosome number to a different column
    df["chrom_num"] = df["chromosome"].apply(lambda x: x[x.find(":") + 1:x.find(":", x.find(":") + 1)])
    # Get the columns to the original order
    cols = df.columns.tolist()
    cols = cols[0:13] + cols[-1:] + cols[13:-1]
    df = df[cols]

    # Filter pseudo-genes
    df = df[df["gene_biotype"] != "polymorphic_pseudogene"]

    # Filter non-coding transcripts
    df = \
    df[df["transcript_biotype"] != "nonsense_mediated_decay"][df["transcript_biotype"] != "non_stop_decay"][
        df["transcript_biotype"] != "polymorphic_pseudogene"]

    df.to_csv(OUTPUT_FILE, sep='\t')
