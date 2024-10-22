import os.path
import pandas as pd
from dsprint.core import COVERAGE_THRESHOLD

INPUT_FILE, COVERAGE_FOLDER = snakemake.input
OUTPUT_FILE = str(snakemake.output)
CHROMOSOME = os.path.splitext(os.path.basename(INPUT_FILE)[len('parsed_chrom'):])[0]


if __name__ == '__main__':

    coverage_df = pd.read_csv(os.path.join(COVERAGE_FOLDER, f'Panel.chr{CHROMOSOME}.coverage.txt.gz'), compression='gzip', sep='\t')
    coverage_df = coverage_df.astype({'#chrom': str})
    coverage_df = coverage_df[coverage_df['#chrom'] == CHROMOSOME]
    coverage_series = coverage_df.set_index('pos')['mean']

    chromosome_df = pd.read_csv(INPUT_FILE, sep='\t', index_col=0)
    chromosome_df = chromosome_df[chromosome_df["FILTER"] == 'PASS']
    merged_df = chromosome_df.merge(coverage_series, how='left', left_on='POS', right_index=True)
    merged_df.rename(columns={'mean': 'COVERAGE'}, inplace=True)

    merged_df = merged_df[merged_df["COVERAGE"] >= COVERAGE_THRESHOLD]
    merged_df.to_csv(OUTPUT_FILE, index=False)
