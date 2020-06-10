import glob
import os.path
import pandas as pd
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 4:
        print('Usage: <script> <input_folder> <chrom_coverage_csv_folder> <output_folder>')
        sys.exit(0)

    INPUT_FOLDER, COVERAGE_FOLDER, OUTPUT_FOLDER = sys.argv[1:]
else:
    INPUT_FOLDER = snakemake.input.input_folder
    COVERAGE_CSVS = snakemake.input.coverage_csvs
    OUTPUT_FOLDER, = snakemake.output


# Note: This script does not operate on a domain-by-domain basis because it is too expensive
# (both in time and memory) to load the chromosome-specific coverage information.
# It thus loads all that coverage data once in memory, and then goes through all domains.

def add_coverage(input_pik, coverage_df_dict, output_folder, skip_existing=True):

    domain = os.path.basename(input_pik).split('.')[0]
    output_file = os.path.join(output_folder, f'{domain}.pik')

    if skip_existing and os.path.exists(output_file):
        print(f'File exists ({output_file}. Skipping..')
        return

    print(f'Processing domain {domain}')
    with open(input_pik, 'rb') as f:
        states_dict = pickle.load(f)

    for state, ds in states_dict.items():
        for d in ds:
            chromosome = d['chrom']
            if chromosome not in coverage_df_dict:
                continue
            df = coverage_df_dict[chromosome]
            df = df[df.index.isin(d['chrom_pos'])]
            d['coverage'] = 0 if df.empty else df.COVERAGE.mean()

    with open(output_file, 'wb') as f:
        pickle.dump(states_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':

    if not os.path.exists(OUTPUT_FOLDER):
        os.mkdir(OUTPUT_FOLDER)

    dfs = {}

    print('Loading coverage data')
    for coverage_csv in COVERAGE_CSVS:
        print(f'  {coverage_csv}')

        df = pd.read_csv(coverage_csv, sep=',', dtype={'CHROM': str}, usecols=['CHROM', 'POS', 'COVERAGE'])
        df = df.set_index('POS')

        chromosome = str(df.iloc[0].CHROM)
        dfs[chromosome] = df.drop('CHROM', axis=1)

    print('Coverage data loaded')

    with ThreadPoolExecutor(max_workers=7) as executor:
        futures = {executor.submit(add_coverage, input_pik, dfs, OUTPUT_FOLDER, True): input_pik
                   for input_pik in glob.glob(f'{INPUT_FOLDER}/*.pik')}
        for future in as_completed(futures):
            print(f'Finished processing {futures[future]}')
