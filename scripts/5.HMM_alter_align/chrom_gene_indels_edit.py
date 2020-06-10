import os
import pandas as pd
from indels_func import table_editing
import pickle
from dsprint.core import INSTANCE_THRESHOLD
from dsprint.mapping_func import create_exon_pos_table


try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 7:
        print('Usage: <script> <domain_states_df> <chromosome_path> <canonic_prot_folder> <hmm_folder> <frameshift_file> <output_path>')
        sys.exit(0)

    DOMAIN_STATES_DF, CHROMOSOME_PATH, CANONIC_PROT_FOLDER, HMM_FOLDER, FRAMESHIFT_FILE, OUTPUT_PATH = sys.argv[1:]
else:
    DOMAIN_STATES_DF, CHROMOSOME_PATH, CANONIC_PROT_FOLDER, HMM_FOLDER, FRAMESHIFT_FILE = snakemake.input
    OUTPUT_PATH = str(snakemake.output[0])

CHROMOSOME = os.path.splitext(os.path.basename(CHROMOSOME_PATH)[len('parsed_filtered_chrom'):])[0]

if __name__ == '__main__':

    domain_stats_df = pd.read_csv(DOMAIN_STATES_DF, sep='\t', index_col=0)
    all_domains_list = domain_stats_df.index.tolist()

    chrom_csv = pd.read_csv(CHROMOSOME_PATH, index_col=0)
    chrom_csv = chrom_csv.sort_values(by=['POS'])
    chrom_csv = chrom_csv.reset_index(drop=True)
    chrom_csv.fillna('', inplace=True)
    chrom_csv['COMMENTS'] = ''

    os.makedirs(OUTPUT_PATH, exist_ok=True)

    for domain_name in all_domains_list:

        with open(os.path.join(CANONIC_PROT_FOLDER, f'{domain_name}_canonic_prot.pik'), 'rb') as f:
            canonic_protein = pickle.load(f)

        domain_data = pd.read_csv(os.path.join(HMM_FOLDER, f'{domain_name}.csv'), sep='\t', index_col=0, dtype={"chrom_num": str})
        domain_data = domain_data[domain_data["chrom_num"] == CHROMOSOME]

        for ens_gene in domain_data["gene"].unique():

            canonic_prot = canonic_protein[ens_gene]
            canonic_prot_t = canonic_prot[:canonic_prot.find(".")]  # Trimming the ".#" at the end
            domain_gene_table = domain_data[domain_data["prot"] == canonic_prot]

            # Making sure that if two HMM-matches overlaps, the higher bit score will come first in the table
            domain_gene_table = domain_gene_table.sort_values(by="BitScore", ascending=False)
            domain_gene_name = domain_gene_table["hugoSymbol"].unique()[0]
            if len(domain_gene_table["hugoSymbol"].unique()) > 1:
                raise RuntimeError(ens_gene + ": more than one Hugo symbol")

            # Creating a table of the exons for this gene, according to the canonical protein
            chrom_raw_data = domain_gene_table["chromosome"].unique()[0]  # there should be only one element here
            if len(domain_gene_table["chromosome"].unique()) > 1:
                raise RuntimeError(ens_gene + ": more than one chromosome raw data")
            targetid = domain_gene_table["#TargetID"].unique()[0]
            exon_table = create_exon_pos_table(chrom_raw_data, targetid, FRAMESHIFT_FILE)

            # Filtering the chromosome data to the gene exons region

            # in case of complement, the minimal position could be at the last row
            exons_start_pos = min(exon_table["start_pos"][0], exon_table["start_pos"][len(exon_table) - 1])

            # in case of complement, the maximal position could be at the first row
            exons_end_pos = max(exon_table["end_pos"][0], exon_table["end_pos"][len(exon_table) - 1])

            chrom_gene_table = \
            chrom_csv[chrom_csv["POS"] >= int(exons_start_pos)][chrom_csv["POS"] <= int(exons_end_pos)][
                chrom_csv["ENSP"] == canonic_prot_t]
            chrom_gene_table = chrom_gene_table.reset_index(drop=True)

            # Adding chrom column to the table
            chrom_gene_table["CHROM"] = CHROMOSOME

            indels_table = table_editing(chrom_gene_table)

            os.makedirs(os.path.join(OUTPUT_PATH, domain_name, ens_gene), exist_ok=True)
            chrom_gene_table.to_csv(os.path.join(OUTPUT_PATH, domain_name, ens_gene,  'chrom_gene_table.csv'))
            indels_table.to_csv(os.path.join(OUTPUT_PATH, domain_name, ens_gene, 'indels_table.csv'))
            exon_table.to_csv(os.path.join(OUTPUT_PATH, domain_name, ens_gene, 'exon_table.csv'))
