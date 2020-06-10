import os.path
import gzip
import pandas as pd

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <tsv_file> <output_folder>')
        sys.exit(0)

    tsv_file, output_folder = sys.argv[1:]
else:
    tsv_file, = snakemake.input
    output_folder = snakemake.output.output_folder


def save_data(output_folder, data, protein):
    pd.DataFrame(data).to_csv(
        os.path.join(output_folder, f'{protein}.jsd.txt'),
        header=['1-Indexed-ProteinPosition', 'ConservationScore', '1-Indexed-ProteinPositionAndResidue'],
        sep='\t',
        index=False
    )


if __name__ == '__main__':

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    header = []
    data = []  # list of lists
    idx_protein = idx_position = idx_score = 0
    last_protein = None

    for i, line in enumerate(gzip.open(tsv_file, 'rt').readlines()):
        line = line.strip()
        if i < 6:
            continue
        elif i == 6:
            header = line.split('\t')
            idx_protein = header.index('#EnsemblProtID_JSD')
            idx_position = header.index('1-Indexed-ProteinPosition')
            idx_score = header.index('ConservationScore')
            idx_residue = header.index('1-Indexed-ProteinPositionAndResidue')
            continue
        else:
            tokens = line.split('\t')
            protein = tokens[idx_protein].split('_')[0]
            if protein != last_protein:
                if last_protein is not None:
                    save_data(output_folder, data, last_protein)
                data = []
                last_protein = protein
                print(f'{protein}')

            data.append([
                tokens[idx_position],
                tokens[idx_score],
                tokens[idx_residue]
            ])

    save_data(output_folder, data, last_protein)
