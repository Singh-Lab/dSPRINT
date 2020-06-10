import gzip
import pickle
from collections import defaultdict
import numpy as np

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 4:
        print('Usage: <script> <pfam_hmm_file> <domains_log_prob_pik> <domains_prob_pik>')
        sys.exit(0)

    INPUT_FILE, DOMAINS_LOG_PROB_FILE, DOMAINS_PROB_FILE = sys.argv[1:]
else:
    INPUT_FILE = snakemake.input[0]
    DOMAINS_LOG_PROB_FILE = snakemake.output[0]
    DOMAINS_PROB_FILE = snakemake.output[1]


if __name__ == '__main__':

    pfam_hmm_file = INPUT_FILE
    _open = gzip.open if pfam_hmm_file.endswith('.gz') else open

    d = {}
    next_state = 0
    domain_name = ''
    log_prob_dict = defaultdict(list)

    for line in _open(pfam_hmm_file, 'rt'):
        parts = line.strip().split()
        first = parts[0]
        if first == 'NAME':
            if domain_name:
                d[domain_name] = log_prob_dict
            domain_name = parts[1]
            next_state = 0
            log_prob_dict = defaultdict(list)

        if first == 'HMM':
            n_states = len(parts) - 1
            next_state = 1

        if next_state > 0 and first == str(next_state):
            log_prob_dict[next_state] = [float(part) for part in parts[1:1 + n_states]]
            next_state += 1

    d[domain_name] = log_prob_dict

    with open(DOMAINS_LOG_PROB_FILE, 'wb') as f:
        pickle.dump(d, f, protocol=pickle.HIGHEST_PROTOCOL)

    for k, v in d.items():
        d[k] = {_k: 1/np.exp(_v) for _k, _v in v.items()}

    with open(DOMAINS_PROB_FILE, 'wb') as f:
        pickle.dump(d, f, protocol=pickle.HIGHEST_PROTOCOL)
