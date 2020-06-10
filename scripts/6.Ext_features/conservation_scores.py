import os.path
import linecache
import pickle


# Modified binary search to find correct fragment file
# Note: a position in a gap with still return the previous fragment — need to check if position is out of bounds later
#
# cons: "phyloP" or "phastCons"
# chrom: chromosome number
# pos: position number
def bin_search(cons, chrom, pos):
    input_path = "conservation_scores/"+cons+"_frags_txt/chr"+str(chrom)+"/"
    with open(input_path+"index.pik",'rb') as handle:
        index = pickle.load(handle)
    if pos < index[0]:
        return(-1)
    return rec_helper(0,len(index)-1,index,pos)


def rec_helper(lo, hi, a, val):
    if hi < lo:
        return -1
    # Indices will never be large, so computing mid in this way is fine
    mid = (hi+lo) / 2
    a_val = a[mid]
    if val >= a_val and mid == len(a)-1:
        return a_val
    elif a_val <= val < a[mid+1]:
        return a_val
    elif val > a_val:
        return rec_helper(mid+1, hi, a, val)
    else:
        return rec_helper(lo, mid-1, a, val)


try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <hmm_states_dict> <conservation_scores_folder>')
        sys.exit(0)

    HMM_STATES_DICT, CONSV_SCORES_FOLDER = sys.argv[1:]
else:
    HMM_STATES_DICT, CONSV_SCORES_FOLDER = snakemake.input


if __name__ == '__main__':

    domain_name = os.path.splitext(os.path.basename(HMM_STATES_DICT))[0]

    with open(HMM_STATES_DICT, 'rb') as f:
        states_dict = pickle.load(f)

    for state, d in states_dict.items():
        chrom = d['chrom']
        for cons in ('phyloP', 'phastCons'):
            frag_path = f"{CONSV_SCORES_FOLDER}/{cons}_frags_txt/chr{chrom}"
            scores = []
            for pos in d['chrom_pos']:
                start = bin_search(cons, chrom, pos)
                score = linecache.getline(f"{frag_path}/{start}.txt", pos - start + 1).strip()

                if score != "":
                    scores.append(float(score))

            d[cons] = scores

        linecache.clearcache()

    with open(HMM_STATES_DICT, 'wb') as f:
        pickle.dump(states_dict, f, protocol=pickle.HIGHEST_PROTOCOL)
