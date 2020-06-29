import pandas as pd
import os.path
from collections import defaultdict
import pickle
import glob

from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

from dsprint.core import retrieve_exon_seq, CHROMOSOMES
from dsprint.mapping_func import create_exon_pos_table


HMMS_FOLDER = snakemake.input.hmm_folder
CANONIC_PROT_FOLDER = snakemake.input.canonic_prot_folder
HG19_FILE = snakemake.input.hg19_file
FRAMESHIFT_FILE = snakemake.input.exon_len_file
OUTPUT_FILE, OUTPUT_TSV_FILE = snakemake.output


if __name__ == '__main__':

    gene_dict = defaultdict(dict)
    all_canonical_proteins = set()

    for domain_file in glob.glob(HMMS_FOLDER + '/*.csv'):

        domain = os.path.splitext(os.path.basename(domain_file))[0]
        print(domain)
        domain_data = pd.read_csv(domain_file, sep='\t', index_col=0, dtype={"chrom_num": str})

        with open(os.path.join(CANONIC_PROT_FOLDER, domain + "_canonic_prot.pik"), 'rb') as f:
            canonic_protein = pickle.load(f)

        for gene, _domain_data in domain_data.groupby('gene'):

            prot_id = canonic_protein[gene]
            if gene in gene_dict and prot_id in gene_dict[gene]:
                continue

            chrom = _domain_data['chrom_num'].unique()[0]
            if chrom not in CHROMOSOMES:
                continue

            chromosome = _domain_data[_domain_data["#TargetID"] == prot_id]['chromosome'].unique()
            if len(chromosome) > 1:
                print(" Error: " + gene + ": more than one chromosome raw data")  # sanity check
            chromosome = chromosome[0]

            exon_table = create_exon_pos_table(chromosome, prot_id, FRAMESHIFT_FILE)
            dna_seq = ''.join([
                retrieve_exon_seq(row["start_pos"], row["end_pos"], chrom, HG19_FILE,
                                  reverse_complement=chromosome.find('complement') >= 0)
                for _, row in exon_table.iterrows()
            ])

            seq = str(Seq(dna_seq, generic_dna).translate())
            seq_len = len(seq)
            gene_dict[gene][prot_id] = seq

            all_canonical_proteins.add((gene, prot_id, seq_len))

    with open(OUTPUT_FILE, 'wb') as f:
        pickle.dump(gene_dict, f, protocol=pickle.HIGHEST_PROTOCOL)

    genes, proteins, lens = tuple(zip(*[(_s[0], _s[1], _s[2]) for _s in all_canonical_proteins]))
    pd.DataFrame({'gene': genes, 'protein': proteins, 'len': lens}).to_csv(OUTPUT_TSV_FILE, sep='\t', index=False, header=False)
