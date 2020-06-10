import fileinput
import pickle
from collections import defaultdict

MATRIX_FILE = '/media/vineetb/t5-vineetb/dsprint/in/BLOSUM62.txt'
OUTPUT_FILE = 'BLOSUM62_dict.pik'

# Saved dictionary is of the form
# {'A': {'A': 4, 'R': -1 ..}, 'R': {'A': -1, ..} .. }

if __name__ == '__main__':

    blosum62_dict = defaultdict(dict)
    for line in open(MATRIX_FILE, 'r').readlines():
        if line[0] == '#':
            continue

        # Reading aa order
        if line[0] == ' ':
            aa_order = line.split()
            continue

        # Reading values to dict
        vals = line.split()
        for i in range(1, len(vals)):
            blosum62_dict[vals[0]][aa_order[i - 1]] = int(vals[i])

    fileinput.close()

    with open(OUTPUT_FILE, 'wb') as handle:
        pickle.dump(blosum62_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
