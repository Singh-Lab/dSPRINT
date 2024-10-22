import os
import os.path
import pickle
import json
from tempfile import NamedTemporaryFile
import shutil
from Bio.Blast.Applications import NcbipsiblastCommandline

import dsprint

DOMAIN_SEQUENCES_DICT = snakemake.input.domain_sequences_dict
PROCESSED_PSSMS_FOLDER = snakemake.input.preprocessed_pssms_folder
OUTPUT_FOLDER = snakemake.output.output_folder

with open(os.path.join(os.path.dirname(dsprint.__file__), '../config.json'), 'r') as f:
    config = json.load(f)
    BLAST_PATH = config['paths']['blast']['bin']
    BLAST_CONFIG = config['blast']
    REMOTE = BLAST_CONFIG['remote']


if __name__ == '__main__':

    os.makedirs(OUTPUT_FOLDER, exist_ok=True)

    with open(DOMAIN_SEQUENCES_DICT, 'rb') as f:
        domain_sequences_dict = pickle.load(f)

    for domain in domain_sequences_dict:

        n_genes = len(domain_sequences_dict[domain])
        for i, (gene, seq) in enumerate(domain_sequences_dict[domain].items(), start=1):

            out_file_name = os.path.join(PROCESSED_PSSMS_FOLDER, f'{gene}.pssm')
            if not os.path.exists(out_file_name):
                print(f'Generating {out_file_name} using BLAST. This will take a while.'
                      'You may want to pre-generate this file beforehand.')
                seq = seq.replace('*', '')

                with NamedTemporaryFile(mode='w') as in_file:
                    in_file.write(seq)
                    in_file.flush()

                    if REMOTE:
                        db = BLAST_CONFIG['default_db']
                        kwargs = {}
                    else:
                        db = config['paths']['blast']['dbs'][BLAST_CONFIG['default_db']]
                        kwargs = {'num_iterations': BLAST_CONFIG['num_iterations'], 'num_threads': BLAST_CONFIG['num_threads']}

                    cline = NcbipsiblastCommandline(
                        cmd=os.path.join(BLAST_PATH, 'bin', 'psiblast'),
                        query=in_file.name,
                        db=db,
                        num_alignments=1,
                        out_ascii_pssm=out_file_name,
                        remote=REMOTE,
                        **kwargs
                    )
                    print(cline)
                    stdout, stderr = cline()

            # Note: It's possible that out_file_name doesn't get generated at all, if there are no hits found.
            # In this case, generate an empty .pssm file so we don't try to regenerate this file the next time this
            # gene/protein is encountered. We'll also have to check for 0 sized files to make this work of course.
            open(out_file_name, 'a').close()  # no-op if file exists, creates an empty file if it doesn't

            # Copy all 'valid' .pssm files (with size>0) to the output folder for downstream processing.
            if os.path.getsize(out_file_name) > 0:
                shutil.copyfile(out_file_name, os.path.join(OUTPUT_FOLDER, f'{gene}.pssm'))

            print(f'Finished gene {i}/{n_genes}')

        print(f'Finished domain {domain}')
