import pandas as pd
import glob
import os.path
import pickle
from Bio import SeqIO

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 5:
        print('Usage: <script> <hmm_folder> <uniprot_fasta_file> <uniprot_idmapping_file> <output_folder>')
        sys.exit(0)

    HMMS_FOLDER, UNIPROT_FASTA_FILE, UNIPROT_IDMAPPING_FILE, OUTPUT_FOLDER = sys.argv[1:]
else:
    HMMS_FOLDER = snakemake.input[0]
    UNIPROT_FASTA_FILE = snakemake.input[1]
    UNIPROT_IDMAPPING_FILE = snakemake.input[2]
    OUTPUT_FOLDER = str(snakemake.output)


def uniprot_canonical_lengths(uniprot_fasta_file, uniprot_idmapping_file):
    # Uniprot fasta headers, as per https://www.uniprot.org/help/fasta-headers, are:
    # >db|UniqueIdentifier|EntryName ProteinName OS=OrganismName ..
    canonical_sequences = {}
    for seq in SeqIO.parse(uniprot_fasta_file, 'fasta'):
        canonical_sequences[seq.name.split('|')[1]] = str(seq.seq)

    df = pd.read_csv(uniprot_idmapping_file, sep='\t', header=None)
    df.columns = ['UniProtKB_AC', 'ID_type', 'ID']
    df = df[df.ID_type == 'Ensembl']

    df['canonical_len'] = df.apply(lambda row: len(canonical_sequences.get(row.UniProtKB_AC, '')), axis=1)
    df = df[df.canonical_len > 0]

    # In some cases, multiple entries exist for a single gene id (representing isoforms)
    # In these cases, take the one representing the longest protein sequence
    df = df.loc[df.groupby('ID', as_index=False)['canonical_len'].idxmax()]
    # We should now have unique IDs, which we can set as the table index.
    df = df.set_index('ID', verify_integrity=True)

    return df


if __name__ == '__main__':

    if not os.path.exists(OUTPUT_FOLDER):
        os.mkdir(OUTPUT_FOLDER)

    UNIPROT = uniprot_canonical_lengths(UNIPROT_FASTA_FILE, UNIPROT_IDMAPPING_FILE)
    domains_files = glob.glob(HMMS_FOLDER + '/*.csv')

    for dom_filename in domains_files:
        domain = os.path.splitext(os.path.basename(dom_filename))[0]
        print(domain)
        domain_data = pd.read_csv(dom_filename, sep='\t', index_col=0)
        canonic_protein = {}

        for gene_id, gene_table in domain_data.groupby('gene'):
            protein_ids = gene_table['prot'].unique()

            if len(protein_ids) == 1:
                canonic_protein[gene_id] = protein_ids[0]

            # If there's more then one transcript: find the canonic protein length from Uniprot tables
            else:
                gene_id_major = gene_id.split('.')[0]
                if gene_id_major not in UNIPROT.index:  # Uniprot ID unavailable? Just take the longest protein
                    canonic_protein[gene_id] = gene_table.loc[gene_table.length.idxmax()].prot
                else:
                    canonic_len = UNIPROT.loc[gene_id_major].canonical_len
                    if canonic_len == 0:  # No sequence information for Uniprot ID? Just take the longest protein
                        canonic_protein[gene_id] = gene_table.loc[gene_table.length.idxmax()].prot
                    else:
                        _gene_table = gene_table[gene_table.length == canonic_len]
                        if not _gene_table.empty:
                            canonic_protein[gene_id] = _gene_table.prot.iloc[0]
                        else:  # No protein matching canonical length? Just take the longest protein
                            canonic_protein[gene_id] = gene_table.loc[gene_table.length.idxmax()].prot

        with open(os.path.join(OUTPUT_FOLDER, domain + '_canonic_prot.pik'), 'wb') as f:
            pickle.dump(canonic_protein, f, protocol=pickle.HIGHEST_PROTOCOL)
