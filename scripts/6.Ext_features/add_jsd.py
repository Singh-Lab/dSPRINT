import sys
import os.path
import glob
import pandas as pd
import pickle
from functools import lru_cache

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 6:
        print('Usage: <script> <hmm_states_folder> <canonic_prot_folder> <jsd_folder> <is_legacy> <output_folder>')
        sys.exit(0)

    HMM_STATES_FOLDER, CANONIC_PROT_FOLDER, JSD_FOLDER, LEGACY, OUTPUT_FOLDER = sys.argv[1:]
    LEGACY = LEGACY == 'True'
else:
    LEGACY = snakemake.params.legacy
    HMM_STATES_FOLDER, CANONIC_PROT_FOLDER, JSD_FOLDER = snakemake.input
    OUTPUT_FOLDER, = snakemake.output


@lru_cache(maxsize=8192)
def fetch_jsd_table(chromosome, gene, protein):
    filepath = f'{JSD_FOLDER}/{chromosome}/{gene}.*/{protein}.*.jsd.txt' if LEGACY else f'{JSD_FOLDER}/{protein}.*.jsd.txt'
    try:
        jsd_file = glob.glob(filepath)[0]
    except IndexError:
        print(f'{protein} (gene={gene}, chrom={chromosome}) file missing from jsd scores')
        return None
    else:
        df = pd.read_csv(jsd_file, header=None, skiprows=1, sep='\t')
        if LEGACY:
            df.columns = ['position', 'aa', 'score']
            df['position'] = df['position'] + 1  # make 1-indexed
        else:
            df.columns = ['position', 'score', 'position_aa']
        return df


if __name__ == '__main__':

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    for pik_file in glob.glob(f'{HMM_STATES_FOLDER}/*.pik'):
        domain = os.path.splitext(os.path.basename(pik_file))[0]

        with open(os.path.join(CANONIC_PROT_FOLDER, f'{domain}_canonic_prot.pik'), 'rb') as f:
            canonic_protein = pickle.load(f)

        with open(pik_file, 'rb') as f:
            states_dict = pickle.load(f)

        for state, ds in states_dict.items():

            for d in ds:
                chromosome = d['chrom']
                gene = d['ens_gene']
                protein_position = d['prot_pos']
                states_dict_aa = d['aa_ref_orig']

                protein = canonic_protein[gene]

                # The jsd gene/protein may have a different version, so strip those off before looking
                gene = gene.split('.')[0]
                protein = protein.split('.')[0]

                if LEGACY:

                    jsd_table = fetch_jsd_table(chromosome, gene, protein)
                    if jsd_table is None:
                        continue

                    # jsd positions are 0-indexed
                    try:
                        row = jsd_table[jsd_table['position'] == protein_position - 1].iloc[0]
                    except IndexError:
                        print(f'{protein} (gene={gene}, chrom={chromosome}) position {protein_position-1} missing')
                        sys.exit(1)
                    jsd_table_aa = row.aa

                else:

                    jsd_table = fetch_jsd_table(chromosome, gene, protein)
                    if jsd_table is None:
                        continue

                    # jsd positions are 1-indexed
                    try:
                        row = jsd_table[jsd_table['position'] == protein_position].iloc[0]
                    except IndexError:
                        print(f'{protein} (gene={gene}, chrom={chromosome}) position {protein_position} missing')
                        sys.exit(1)
                    jsd_table_aa = row.position_aa.split(':')[1]

                if jsd_table_aa != "-" and jsd_table_aa.upper() != states_dict_aa.upper():
                    if not (
                            jsd_table_aa == '*' and states_dict_aa.upper() == 'X'):  # mismatch isn't stop codon notation
                        d['100-way-BLOSUM_JSD'] = -1
                else:
                    d['100-way-BLOSUM_JSD'] = row.score

        with open(os.path.join(OUTPUT_FOLDER, f'{domain}.pik'), 'wb') as f:
            pickle.dump(states_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
