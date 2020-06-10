import os.path
import pickle
from collections import defaultdict

from dsprint.core import CHROMOSOMES

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) !=3 :
        print('Usage: <script> <exon_shift_folder> <output_file>')
        sys.exit(0)

    INPUT_FOLDER, OUTPUT_FILE = sys.argv[1:]
else:
    INPUT_FOLDER = snakemake.input[0]
    OUTPUT_FILE = snakemake.output[0]


if __name__ == '__main__':

    d = defaultdict(list)

    for chromosome in CHROMOSOMES:
        for gene in os.listdir(os.path.join(INPUT_FOLDER, chromosome)):
            gene_dir = os.path.join(INPUT_FOLDER, chromosome, gene)
            for exon_file in os.listdir(gene_dir):
                for line in open(os.path.join(gene_dir, exon_file), 'r'):
                    line = line.strip()
                    if not line.startswith('>'):
                        parts = line.split()
                        start, end = map(int, parts[0].split(':'))
                        if start == end and start < 0:
                            d[exon_file].append((-start, len(parts[1]), parts[1]))

    with open(OUTPUT_FILE, 'wb') as f:
        pickle.dump(d, f, protocol=pickle.HIGHEST_PROTOCOL)
