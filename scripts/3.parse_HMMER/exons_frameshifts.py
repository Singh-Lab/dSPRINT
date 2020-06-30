import os.path
from glob import glob
import pickle
from collections import defaultdict

INPUT_FOLDERS = snakemake.input
OUTPUT_FILE = snakemake.output[0]


if __name__ == '__main__':

    d = defaultdict(list)

    for input_folder in INPUT_FOLDERS:
        print(f'Examining genes in folder {input_folder}')
        for gene in os.listdir(input_folder):
            gene_dir = os.path.join(input_folder, gene)
            if not os.path.isdir(gene_dir):
                continue
            for exon_file in glob(f'{gene_dir}/*.exons.txt'):
                for line in open(exon_file, 'r'):
                    line = line.strip()
                    if not line.startswith('>'):
                        parts = line.split()
                        start, end = map(int, parts[0].split(':'))
                        if start == end and start < 0:
                            d[exon_file].append((-start, len(parts[1]), parts[1]))

    with open(OUTPUT_FILE, 'wb') as f:
        pickle.dump(d, f, protocol=pickle.HIGHEST_PROTOCOL)
