import os.path
import pickle
import pandas as pd

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) < 3:
        print('Usage: <script> <pfam_clan_csv> <output_folder> [<human_proteome_csv>]')
        sys.exit(0)

    INPUT_FILE, OUTPUT_FOLDER = sys.argv[1:3]
    PROTEOME_CSV = sys.argv[3] if len(sys.argv) == 4 else None
else:
    INPUT_FILE = snakemake.input[0]
    PROTEOME_CSV = snakemake.input[1]
    OUTPUT_FOLDER = os.path.dirname(snakemake.output[0])


if __name__ == '__main__':
    df = pd.read_csv(INPUT_FILE, sep='\t', header=None)
    df.columns = ['pfam_id', 'clan_id', 'clan_name', 'domain_name', 'description']

    # Each domain_name belongs to a single clan_name
    with open(os.path.join(OUTPUT_FOLDER, 'updated_domain_to_clan_dict.pik'), 'wb') as f:
        pickle.dump(
            df.set_index('domain_name', verify_integrity=True)['clan_name'].to_dict(),
            f, protocol=pickle.HIGHEST_PROTOCOL
        )

    # Each clan_name can have many domain_names
    with open(os.path.join(OUTPUT_FOLDER, 'updated_clan_to_domains_dict.pik'), 'wb') as f:
        pickle.dump(
            df.set_index('clan_name')['domain_name'].groupby(level=0).agg(list).to_dict(),
            f, protocol=pickle.HIGHEST_PROTOCOL
        )

    # Each domain_name has a single pfam_id
    with open(os.path.join(OUTPUT_FOLDER, 'updated_domain_to_pfam_acc_dict.pik'), 'wb') as f:
        pickle.dump(
            df.set_index('domain_name', verify_integrity=True)['pfam_id'].to_dict(),
            f, protocol=pickle.HIGHEST_PROTOCOL
        )

    # Extra data for human proteome
    # This section generates:
    #    domain_to_clan_dict.pik
    #    clan_to_domains_dict.pik
    #    domain_to_pfam_acc_dict.pik

    if PROTEOME_CSV is not None:
        proteome_df = pd.read_csv(PROTEOME_CSV, sep='\t', skiprows=[0, 1, 2], index_col=False, header=None)
        proteome_header = pd.read_csv(PROTEOME_CSV, sep='<', skiprows=[0, 1], nrows=1, header=None)
        proteome_header = proteome_header.iloc[0].tolist()
        proteome_header.remove("#")
        proteome_header = [x[:x.find('>')] for x in proteome_header]
        proteome_df.columns = proteome_header

        # Note: 'hmm name' here corresponds to the domain name
        # Each domain name belongs to a single clan_name
        with open(os.path.join(OUTPUT_FOLDER, 'domain_to_clan_dict.pik'), 'wb') as f:
            # Note: Each 'hmm name' can appear multiple times in the csv
            pickle.dump(
                proteome_df.set_index('hmm name', verify_integrity=False)['clan'].to_dict(),
                f, protocol=pickle.HIGHEST_PROTOCOL
            )

        # Each clan_name can have many domain names
        with open(os.path.join(OUTPUT_FOLDER, 'clan_to_domains_dict.pik'), 'wb') as f:
            pickle.dump(
                proteome_df.set_index('clan')['hmm name'].groupby(level=0).agg(list).to_dict(),
                f, protocol=pickle.HIGHEST_PROTOCOL
            )

        # Each domain name has an accession string
        with open(os.path.join(OUTPUT_FOLDER, 'domain_to_pfam_acc_dict.pik'), 'wb') as f:
            # Note: Each 'hmm name' can appear multiple times in the csv
            pickle.dump(
                proteome_df.set_index('hmm name', verify_integrity=False)['clan'].to_dict(),
                f, protocol=pickle.HIGHEST_PROTOCOL
            )
