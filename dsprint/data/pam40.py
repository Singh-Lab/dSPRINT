from collections import defaultdict
import fileinput
import pickle

MATRIX_FILE = '/media/vineetb/t5-vineetb/dsprint/in/pam40.txt'
OUTPUT_FILE = 'PAM40_dict.pik'

if __name__ == '__main__':

    pam40_dict = defaultdict(dict)
    for line in open(MATRIX_FILE, 'r').readlines():
        # Reading aa order
        if line[0] == " ":
            aa_order = line.split()
            continue

        # Reading values to dict
        vals = line.split()
        for i in range(1, len(vals)):
            pam40_dict[vals[0]][aa_order[i - 1]] = int(vals[i])

    fileinput.close()

    # Saving to file
    with open(OUTPUT_FILE, 'wb') as handle:
        pickle.dump(pam40_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)