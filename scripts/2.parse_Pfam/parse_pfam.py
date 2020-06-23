import gzip
import pandas as pd

# Fields that we are interested in - for multiple valued fields (e.g. GA), the first component is extracted
FIELDS = 'NAME', 'LENG', 'GA'

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


