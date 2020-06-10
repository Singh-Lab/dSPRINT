import gzip
import pandas as pd

# Fields that we are interested in - for multiple valued fields (e.g. GA), the first component is extracted
FIELDS = 'NAME', 'LENG', 'GA'

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <pfam_hmm_file> <output_csv_file>')
        sys.exit(0)

    INPUT_FILE, OUTPUT_FILE = sys.argv[1:]
else:
    INPUT_FILE = snakemake.input[0]
    OUTPUT_FILE = snakemake.output[0]


if __name__ == '__main__':

    pfam_hmm_file = INPUT_FILE
    _open = gzip.open if pfam_hmm_file.endswith('.gz') else open

    d = {k: [] for k in FIELDS}
    for line in _open(pfam_hmm_file, 'rt'):
        parts = line.strip().split()
        first = parts[0]
        if first in FIELDS:
            d[first].append(parts[1])

    pd.DataFrame(d).set_index('NAME').to_csv(OUTPUT_FILE, sep='\t')


