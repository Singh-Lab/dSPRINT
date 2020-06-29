import sys
import os.path
from types import SimpleNamespace
import runpy
from dsprint.core import CHROMOSOMES


SCRIPT_PATH = 'scripts/6.Ext_features/add_spider2.py'

if __name__ == '__main__':
    snakemake = SimpleNamespace(
        input='_out/hmms,_out/hmm_states_1,_out/domains_canonic_prot,_out/spd3,_out/all_domains_genes_prot_seq.pik'.split(','),
        output=['_out/hmm_states_2'],
        params=SimpleNamespace(legacy=False)
    )

    sys.path.append(os.path.dirname(SCRIPT_PATH))
    runpy.run_path(path_name=SCRIPT_PATH, init_globals={'snakemake': snakemake}, run_name='__main__')
