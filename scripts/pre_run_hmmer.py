"""
run-hmmer insists that all individual domains that it processes be named in a very specific way, and located in very
specific folders. This script prepares such for run-hmmer's consumption.

Note that though run-hmmer is capable of creating these files, it only does so when downloading new PFAM data.
"""
import os.path

INPUT_FILE = snakemake.input[0]
OUTPUT_FOLDER = snakemake.output[0]


if __name__ == '__main__':

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    lines = []
    for l in open(INPUT_FILE, 'r'):
        lines.append(l)

        if l.strip() == '//':
            if len(lines) > 0:
                f = open(os.path.join(OUTPUT_FOLDER, f'{accession}_{name}.hmm'), 'w')
                f.writelines(lines)
                f.close()
            name = accession = ''
            lines = []

        else:
            if l.startswith('NAME '):
                name = l.strip().split()[-1]
            elif l.startswith('ACC '):
                accession = l.strip().split()[-1].split('.')[0]
