import pandas as pd
import numpy as np


INPUT_CSV = snakemake.input.input_csv
OUTPUT_CSV = snakemake.output.output_csv


def calc_relevant_idx(pos_idx, window_size, max_idx):
    """Calculate the relevant domain positions for the window size"""

    idx_list = [pos_idx]

    for i in range(1, window_size + 1):

        if (pos_idx - i) > 0:
            idx_list.append(pos_idx - i)
        if (pos_idx + i) <= max_idx:
            idx_list.append(pos_idx + i)

    idx_list.sort()
    return idx_list


def calc_windowed_feature(domain_name, feature_name, window_size, domains_features_df):
    """Calculate a windowed (taking the mean across the window) feature for the input domain.
    Returning a list of the windowed feature in the order of the domain positions in the input table."""

    curr_domain_table = domains_features_df[domains_features_df["domain_name"] == domain_name]
    max_pos = max([int(index[index.rfind("_") + 1:]) for index in curr_domain_table.index.tolist()])

    # init features_lists
    feature_domain_mean_values = []
    feature_domain_std_values = []

    for index, row in curr_domain_table.iterrows():
        curr_pos = int(index[index.rfind("_") + 1:])
        window_idx = calc_relevant_idx(curr_pos, window_size, max_pos)

        # Add relevant feature values from the positions in the window
        curr_pos_feature_list = []
        for pos in window_idx:
            idx = domain_name + "_" + str(pos)
            try:
                feature_val = curr_domain_table.loc[idx, :][feature_name]
                curr_pos_feature_list.append(feature_val)
            except:
                # The relevant idx isn't present, don't add it's value
                continue

        feature_domain_mean_values.append(np.mean(curr_pos_feature_list))
        feature_domain_std_values.append(np.std(curr_pos_feature_list))

    return feature_domain_mean_values, feature_domain_std_values, max_pos


if __name__ == '__main__':

    domains_features_df = pd.read_csv(INPUT_CSV, sep='\t', index_col=0)
    domains_list = domains_features_df["domain_name"].unique().tolist()

    WINDOWS = [1, 3, 5, 10]
    features = (
        "avg_maf_altered",
        "phyloP1_avg",
        "phyloP2_avg",
        "phyloP3_avg",
        "blosum_avg",
        "pam_avg",
        "pseudo_dNdS",
        "pfam_prob_max",
        "sift_avg",
        "polyphen_avg",
        "avg_clinvar_score",
        "med_jsd_100way_blosum",
        "jsds_ratio",
        "hindex_avg",
        "vol_avg",
        "aa_ref_charge_majority",
        "aa_ref_alpha_prop_avg",
        "aa_ref_beta_prop_avg",
        "aa_ref_turn_prop_avg" ,
        "H_bond_donor_avg",
        "H_bond_acceptor_avg",
        "sub_diff_hindex_avg_weighted",
        "sub_diff_vol_avg_weighted",
        "sub_func_group_stay_freq",
        "sub_func_group_move_freq",
        "solvent_acc_avg",
        "solvent_acc_std",
        "hsa2_cn_avg",
        "hsb2_cn_avg",
        "backbone_Phi_angle_avg",
        "backbone_Psi_angle_avg",
        "c-alpha_tau_angle_avg",
        "c-alpha_theta_angle_avg",
        "helix_prob_avg",
        "sheet_prob_avg",
        "turn_prob_avg",
        "hsa2_HSE-up_avg",
        "hsa2_HSE-down_avg",
        "hsb2_HSE-up_avg",
        "hsb2_HSE-down_avg"
    )

    for feature in features:
        for window_size in WINDOWS:
            feature_vals_mean_list = []
            feature_vals_std_list = []
            len_vals_list = []

            for domain_name in domains_list:
                domain_feature_mean_vals, domain_feature_std_vals, _ = calc_windowed_feature(domain_name,
                                                                                                     feature,
                                                                                                     window_size,
                                                                                                     domains_features_df)
                feature_vals_mean_list.extend(domain_feature_mean_vals)
                feature_vals_std_list.extend(domain_feature_std_vals)

            domains_features_df[f'wm-{window_size}-{feature}'] = feature_vals_mean_list
            domains_features_df[f'ws-{window_size}-{feature}'] = feature_vals_std_list

    domains_features_df.to_csv(OUTPUT_CSV, sep='\t')
