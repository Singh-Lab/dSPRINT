import pandas as pd
import pickle
import os.path


pfam_aa_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

try:
    snakemake
except NameError:
    import sys

    if len(sys.argv) != 5:
        print('Usage: <script> <hmmer_results> <pfam_data> <prob_dict> <output_folder>')
        sys.exit(0)

    HMMER_RESULTS, PFAM_DATA, PROB_DICT, OUTPUT_FOLDER = sys.argv[1:]
else:
    HMMER_RESULTS = snakemake.input[0]
    PFAM_DATA = snakemake.input[1]
    PROB_DICT = snakemake.input[2]
    OUTPUT_FOLDER = str(snakemake.output)


def domain_conserved_states_filter(domain_data, domain_hmm_prob, con_threshold):
    """
    Filter the given domain data to the domain instances that contain the major allele of the conserved states.
    Conserved states are determined as having Pfam hmm emission prob. above the "con_threshold" given.
    """
    # Find the conserved states and their major allele:
    con_states_dict = {}
    for state in domain_hmm_prob.keys():
        prob_list = domain_hmm_prob[state]
        for i in range(len(prob_list)):
            p = prob_list[i]
            if p >= con_threshold:
                major_allele = pfam_aa_order[i]
                con_states_dict[state] = major_allele

    # Filter the domain instances missing any conserved state
    # Creating a new data frame for the filtered table
    if len(con_states_dict.keys()) > 0:

        # Creating a new data frame for the filtered table
        domain_filtered = pd.DataFrame(columns=domain_data.columns)
        domain_filtered_i = 0

        # Iterating over all domain instances and check the conserved states
        for index, row in domain_data.iterrows():
            target_seq = list(row['Target_Seq'])
            hmm_pos = (row['HMM_Pos']).split(',')
            add_flag = True
            for con_state in con_states_dict.keys():
                try:
                    hmm_pos.index(str(con_state))
                # The conserved state is not even in the alignment
                except ValueError:
                    add_flag = False
                    break
                state_idx = hmm_pos.index(str(con_state))
                aligned_aa = target_seq[state_idx]

                # Compare the instance aa with the major allele
                if aligned_aa != con_states_dict[con_state]:
                    add_flag = False
                    break

            if add_flag:
                new_row = row.copy(deep=True)
                domain_filtered.loc[domain_filtered_i] = new_row
                domain_filtered_i += 1

        return domain_filtered
    else:
        return domain_data


if __name__ == '__main__':

    if not os.path.exists(OUTPUT_FOLDER):
        os.mkdir(OUTPUT_FOLDER)

    allhmm = pd.read_csv(HMMER_RESULTS, sep='\t', index_col=0)
    pfam_df = pd.read_csv(PFAM_DATA, sep='\t', index_col=0)
    with open(PROB_DICT, 'rb') as f:
        domains_hmm_prob_dict = pickle.load(f)

    # Conservation threshold (after witnessing a small bump in the prob. distribution after 0.99)
    con_threshold = 0.99

    for domain, _df in allhmm.groupby(['domain_name']):

        if domain not in pfam_df.index:
            continue

        domain_details = pfam_df.loc[domain]
        _df = _df[_df['BitScore'] >= domain_details.GA]
        _df = _df[(_df['HMMStart'] == 1) & (_df['HMMEnd'] == int(domain_details.LENG))]

        domain_hmm_prob = domains_hmm_prob_dict[domain]
        _df = domain_conserved_states_filter(_df, domain_hmm_prob, con_threshold)

        if not _df.empty:
            _df.to_csv(os.path.join(OUTPUT_FOLDER, domain + '.csv'), sep='\t')
