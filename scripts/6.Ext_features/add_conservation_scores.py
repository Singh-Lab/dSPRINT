import os
import os.path
import glob
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

from dsprint.conservation_score import ConservationScore

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 6:
        print('Usage: <script> <input_pik_folder> <score_gz_folder> <index_folder> <cons_name> <output_pik_folder>')
        sys.exit(0)

    INPUT_PIK_FOLDER, GZ_FOLDER, INDEX_FOLDER, CONSERVATION_NAME, OUTPUT_PIK_FOLDER = sys.argv[1:]
else:
    CONSERVATION_NAME = snakemake.params.conservation_name
    CHROMOSOMES = snakemake.params.chromosomes
    WIGFIX_FILES = snakemake.input.wigfix_files
    INDEX_FILES = snakemake.input.index_files
    INPUT_PIK_FOLDER = snakemake.input.input_pik_folder
    OUTPUT_PIK_FOLDER, = snakemake.output


def add_conservation_score(input_pik, conservation, output_folder):
    domain = os.path.basename(input_pik).split('.')[0]
    with open(input_pik, 'rb') as f:
        states_dict = pickle.load(f)

    all_scores = []
    for state, ds in states_dict.items():
        # underscored keys indicate an attribute applicable to the whole domain, not a numerical 'state'
        if str(state).startswith('_'):
            continue
        for d in ds:
            chromosome = d['chrom']
            scores = []
            for position in d['chrom_pos']:
                try:
                    score = conservation[chromosome][position]
                except (KeyError, IndexError):
                    pass
                else:
                    scores.append(score)
                    all_scores.append(score)
            d[conservation.name] = scores

    # Note: Important to promote dtypes before calculating mean/std
    # see https://numpy.org/doc/stable/reference/generated/numpy.mean.html etc.
    states_dict[f'_{conservation.name}_mean'] = np.mean(all_scores, dtype=np.float64)
    states_dict[f'_{conservation.name}_std'] = np.std(all_scores, dtype=np.float64)

    with open(os.path.join(output_folder, f'{domain}.pik'.format(domain=domain)), 'wb') as f:
        pickle.dump(states_dict, f, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':

    os.makedirs(OUTPUT_PIK_FOLDER, exist_ok=True)

    conservation_score = ConservationScore(name=CONSERVATION_NAME)
    conservation_score.load(CHROMOSOMES, WIGFIX_FILES, INDEX_FILES)

    with ThreadPoolExecutor(max_workers=1) as executor:
        futures = {
            executor.submit(add_conservation_score, input_pik, conservation_score, OUTPUT_PIK_FOLDER):
                input_pik
            for input_pik in glob.glob(f'{INPUT_PIK_FOLDER}/*.pik')
        }
        for future in as_completed(futures):
            print(f'Finished processing {futures[future]}')
