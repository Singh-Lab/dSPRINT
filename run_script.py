import sys
import os.path
from types import SimpleNamespace
import runpy
from dsprint.core import CHROMOSOMES


SCRIPT_PATH = 'scripts/6.Ext_features/process_jsd_data.py'

if __name__ == '__main__':
    snakemake = SimpleNamespace(
        input=['/media/vineetb/t5-vineetb/dsprint/in/processed/pertinint/100way-jsdconservation_domainweights-GRCh37.txt.gz'],
        output=['_out/jsd_scores']
    )

    sys.path.append(os.path.dirname(SCRIPT_PATH))
    runpy.run_path(path_name=SCRIPT_PATH, init_globals={'snakemake': snakemake}, run_name='__main__')
