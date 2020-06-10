import os
import pandas as pd
import pickle
import glob


def count_domain_instances(domain_gene_table, count_overlaps=True):
    if count_overlaps:
        return domain_gene_table.shape[0]

    else:
        instance_counter = 0
        last_target_end = 0

        for i, row in domain_gene_table.iterrows():
            curr_target_start = int(row["TargetStart"])
            curr_target_end = int(row["TargetEnd"])
            if curr_target_start > last_target_end:
                instance_counter += 1
                last_target_end = curr_target_end
            # If the instance overlpas the previous one
            else:
                # Updating to the smaller traget end
                if curr_target_end > last_target_end:
                    last_target_end = curr_target_end
                # Continue without incrememnting the counter
                continue

        return instance_counter


try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 4:
        print('Usage: <script> <hmms_folder> <canonic_prot_folder> <output_domain_stats_df>')
        sys.exit(0)

    HMMS_FOLDER, CANONIC_PROT_FOLDER, OUTPUT_FILE = sys.argv[1:]
else:
    HMMS_FOLDER = snakemake.input.hmm_folder
    CANONIC_PROT_FOLDER = snakemake.input.canonic_prot_folder
    OUTPUT_FILE = str(snakemake.output[0])


if __name__ == '__main__':

    domains_stats = {}

    for dom_filename in glob.glob(HMMS_FOLDER + '/*.csv'):

        curr_domain_stats = []

        domain_sym = os.path.splitext(os.path.basename(dom_filename))[0]
        domain_data = pd.read_csv(dom_filename, sep='\t', index_col=0)
        domain_ens_genes = (domain_data["gene"]).unique()

        with open(os.path.join(CANONIC_PROT_FOLDER, domain_sym + '_canonic_prot.pik'), 'rb') as handle:
            canonic_protein = pickle.load(handle)

        domain_instance_num = 0
        for ens_gene in domain_ens_genes:
            # Filtering the domain data for this gene according to the canonical protein id
            canonic_prot = canonic_protein[ens_gene]
            canonic_prot_t = canonic_prot[:canonic_prot.find(".")]  # Trimming the ".#" at the end
            domain_gene_table = domain_data[domain_data["prot"] == canonic_prot]

            # Count the number of domain instances in this gene
            domain_gene_table = domain_gene_table.sort_values(by=["TargetStart", "BitScore"], ascending=[True, False])
            domain_instance_num += count_domain_instances(domain_gene_table, count_overlaps=True)

        curr_domain_stats.append(len(domain_ens_genes))  # No. of genes
        curr_domain_stats.append(domain_instance_num)    # No. of domain instances

        domains_stats[domain_sym] = curr_domain_stats

    domains_stats_df = pd.DataFrame.from_dict(domains_stats, orient='index')
    domains_stats_df.columns = ["genes", "instances"]
    domains_stats_df = domains_stats_df.sort_values(by=["instances", "genes"], ascending=[False, False])
    domains_stats_df.to_csv(OUTPUT_FILE, sep='\t')

