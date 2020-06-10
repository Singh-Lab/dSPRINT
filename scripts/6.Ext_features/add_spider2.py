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
    if len(sys.argv) != 7:
        print('Usage: <script> <hmm_folder> <hmm_states_folder> <canonic_prot_folder> <spider_folder> <all_domains_genes_prot_seq_file> <output_folder>')
        sys.exit(0)

    HMM_FOLDER, HMM_STATES_FOLDER, CANONIC_PROT_FOLDER, SPIDER_FOLDER, ALL_DOMAINS_GENE_PROT_SEQ_FILE, OUTPUT_FOLDER = sys.argv[1:]
else:
    HMM_FOLDER, HMM_STATES_FOLDER, CANONIC_PROT_FOLDER, SPIDER_FOLDER, ALL_DOMAINS_GENE_PROT_SEQ_FILE = snakemake.input
    OUTPUT_FOLDER, = snakemake.output


# file extensions by spider that have csv data
SPIDER_EXTS = 'spd3', 'hsa2', 'hsb2'
# columns in the above found files that give us protein positions, and can thus serve as DataFrame indices
SPIDER_INDEX_COLS = '#', '#index', '#index'


@lru_cache(maxsize=8192)
def spider_dataframes(gene):
    files = [os.path.join(SPIDER_FOLDER, f'{gene}.{ext}') for ext in SPIDER_EXTS]
    dfs = [None] * len(files)
    for i, file in enumerate(files):
        if os.path.exists(file):
            dfs[i] = pd.read_csv(file, sep='\t', index_col=SPIDER_INDEX_COLS[i])

    assert all(dfs[0].AA == dfs[1].AA) and all(dfs[1].AA == dfs[2].AA), "AAs in merged DFs don't match!"

    merged = dfs[0].merge(
        dfs[1], left_index=True, right_index=True
    ).merge(
        dfs[2], left_index=True, right_index=True, suffixes=('_hsa2', '_hsb2')
    )

    return merged


if __name__ == '__main__':

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    for pik_file in glob.glob(f'{HMM_STATES_FOLDER}/*.pik'):
        domain = os.path.splitext(os.path.basename(pik_file))[0]

        with open(os.path.join(CANONIC_PROT_FOLDER, f'{domain}_canonic_prot.pik'), 'rb') as f:
            canonic_protein = pickle.load(f)

        with open(ALL_DOMAINS_GENE_PROT_SEQ_FILE, 'rb') as f:
            seqs = pickle.load(f)

        with open(pik_file, 'rb') as f:
            states_dict = pickle.load(f)

        for state, ds in states_dict.items():

            for d in ds:
                chromosome = d['chrom']
                gene = d['ens_gene']
                protein_position = d['prot_pos']
                states_dict_aa = d['aa_ref_orig']

                protein = canonic_protein[gene]

                seq = seqs[gene][protein]
                seq_ = seq[0:d['prot_pos']]
                skip = seq_.count('X') + seq_.count('*') + seq_.count('.') + seq_.count('-')
                spider_pos = d['prot_pos'] - skip

                spider_df = spider_dataframes(gene)

                if spider_pos in spider_df.index:
                    row = spider_df.loc[spider_pos]

                    d["spider2-2nd_struct"] = row["SS"]  # secondary structure prediction
                    d["spider2-helix_prob"] = row["P(H)"]  # alpha-Helix prob.
                    d["spider2-sheet_prob"] = row["P(E)"]  # beta-sheet prob.
                    d["spider2-turn_prob"] = row["P(C)"]  # turn prob.
                    d["spider2-angle_Phi"] = row["Phi"]  # backbone_torsion angle
                    d["spider2-angle_Psi"] = row["Psi"]  # backbone_torsion angle
                    d["spider2-angle_tau"] = row["Tau(i-2=>i+1)"]  # c-alpha angle (i-2=>i+1)
                    d["spider2-angle_theta"] = row["Theta(i-1=>i+1)"]  # c-alpha angle (i-1=>i+1)
                    d["spider2-ASA"] = row["ASA"]  # Accessible Surface Area (solvent accessibility)
                    d["spider2-hsa2_HSEu"] = row["HSEu_hsa2"]  # half-sphere exposure Cα-Cα vectors (HSEα-up)
                    d["spider2-hsa2_HSEd"] = row["HSEd_hsa2"]  # half-sphere exposure Cα-Cα vectors (HSEα-down)
                    d["spider2-hsb2_HSEu"] = row["HSEu_hsb2"]  # half-sphere exposure Cα-Cβ vectors (HSEβ-up)
                    d["spider2-hsb2_HSEd"] = row["HSEd_hsb2"]  # half-sphere exposure Cα-Cβ vectors (HSEβ-down)
                    d["spider2-hsa2_CN"] = row["CN_hsa2"]  # contact number for Cα-Cα
                    d["spider2-hsb2_CN"] = row["CN_hsb2"]  # contact number for Cα-Cβ

        with open(os.path.join(OUTPUT_FOLDER, f'{domain}.pik'), 'wb') as f:
            pickle.dump(states_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
