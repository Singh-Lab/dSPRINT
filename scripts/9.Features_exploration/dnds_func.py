import os.path
import pickle
import dsprint.data as data


with open(os.path.join(os.path.dirname(data.__file__), 'codon_ns_table.pik'), 'rb') as f:
    codon_ns_table = pickle.load(f)


def seq_ns(sequence):
    """Given a sequence of nucletides that comprise full codons triplets
    calculate and return N = the total number of nonsynnonymous sites,
    and S = the total number of synnonymous sites
    """

    N = S = 0

    for i in range(0, len(sequence), 3):
        codon = sequence[i:i + 3]
        N += codon_ns_table[codon]["N"]
        S += codon_ns_table[codon]["S"]

    return N, S
