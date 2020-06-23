import pandas as pd
import numpy as np
import pickle
from collections import defaultdict
import os.path
import glob
from scipy.stats import entropy
from Bio.SeqUtils.ProtParamData import kd
from dnds_func import seq_ns
from aa_chemical_properties import aa_charge, aa_charge_dict, aa_functional_group, aa_functional_group_dict, aa_propensity,\
                                    propensity_chou_fasman, aa_volume_group, aa_volume, aa_volume_group_dict, aa_h_bond_donor, aa_h_bond_acceptor
from ext_predictors_codes import sift_codes, polyphen_codes, clinvar_codes
from calc_exac_freq_func import codon_table
from entropy_func import SE_hist, JSD_background, JSD_hist

from dsprint.core import POPULATIONS_ANS, POPULATIONS_ACS
import dsprint.data as data


SIFT_THRESHOLD = 0.05

# Rare SNP thresholds
MAFT_5 = 0.005
MAFT_05 = 0.0005
MAFT_005 = 0.00005

pfam_aa_order = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
AMINO_ACIDS = pfam_aa_order + ['*']

HMM_STATES_FOLDER = snakemake.input.hmm_states_folder
PROB_DICT = snakemake.input.prob_dict
OUTPUT_CSV = snakemake.output.output_csv


def ExAC_MAF_features(sites_aa_num, sites_aa_alter_num, maf_list):

    d = {}

    # avg MAF
    d['avg_maf_all'] = 0 if sites_aa_num == 0 else np.sum(maf_list) / float(sites_aa_num)

    # avg MAF of all the altered sites
    d['avg_maf_altered'] = 0 if sites_aa_alter_num == 0 else np.sum(maf_list) / float(sites_aa_alter_num)

    bins = [0, 0.001, 0.005, 0.01, 0.02, 0.04, 0.06, 0.08, 0.1, 0.2, 0.5]
    non_zero_maf_lst = np.array(maf_list)[np.nonzero(maf_list)[0].tolist()]
    maf_hist = np.histogram(non_zero_maf_lst, bins)[0]

    for i, maf_hist_value in enumerate(maf_hist):
        d['maf_hist_' + str(bins[i]) + '-' + str(bins[i + 1])] = maf_hist_value

    return d


def ExAC_population_features(pop_maf_list, pop_maf_syn_list, pop_maf_nonsyn_list):

    d = {}

    for i in range(len(an_str)):
        # populations total maf avg
        d['maf_' + an_str[i][3:]] = 0 if len(pop_maf_list[i]) == 0 else np.average(pop_maf_list[i])

    for i in range(len(an_str)):
        # populations syn maf avg
        d['maf_syn_' + an_str[i][3:]] = 0 if len(pop_maf_syn_list[i]) == 0 else np.average(pop_maf_syn_list[i])

    for i in range(len(an_str)):
        # populations non-syn maf avg
        d['maf_nonsyn_' + an_str[i][3:]] = 0 if len(pop_maf_nonsyn_list[i]) == 0 else np.average(pop_maf_nonsyn_list[i])

    return d


def ExAC_count_features(sites_aa_num, sites_aa_alter_num, sites_snp_num, sites_snp_alter_num):

    d = {}

    # Feature: number of alterations - aa level (raw and normalized by total number of matched positions)
    d['alter_num_aa'] = sites_aa_alter_num
    d['alter_num_aa_norm'] = 0 if sites_aa_num == 0 else sites_aa_alter_num / float(sites_aa_num)

    # Feature: number of alterations - DNA level (raw and normalized by total number of matched positions)
    d['alter_num_snp'] = sites_snp_alter_num
    d['alter_num_snp_norm'] = 0 if sites_snp_num == 0 else sites_snp_alter_num / float(sites_snp_num)

    # Feature: average number of poymorphisms at one site
    d['avg_aa_polymorphisms'] = 0 if sites_aa_alter_num == 0 else sites_poly_aa_num / float(sites_aa_alter_num)

    # Feature: fraction of altered sites with more than 1 polymorphism
    d['frac_poly_aa'] = 1 if sites_aa_alter_num == 0 else sites_poly_aa_several / float(sites_aa_alter_num)

    return d


def ExAC_rareSNP_features(sites_snp_alter_num, rare_5_num, rare_05_num, rare_005_num):

    # Feature: fraction of rare SNPs (0.5%, 0.05%, 0.005%)
    return {
        'rare_poly_0.5': 0 if sites_snp_alter_num == 0 else rare_5_num / float(sites_snp_alter_num),
        'rare_poly_0.05': 0 if sites_snp_alter_num == 0 else rare_05_num / float(sites_snp_alter_num),
        'rare_poly_0.005': 0 if sites_snp_alter_num == 0 else rare_005_num / float(sites_snp_alter_num)
    }


def conservation_features(phastCons_dict, phyloP_dict):

    d = {}
    positions = 1, 2, 3  # codon positions

    phastCons = np.vstack(phastCons_dict[p] for p in positions)
    phyloP = np.vstack(phyloP_dict[p] for p in positions)

    # conservation scores avg for each codon position
    phastCons_mean = np.nanmean(phastCons, axis=1)
    for p in positions:
        d[f'phastCons{p}_avg'] = phastCons_mean[p-1]

    phyloP_mean = np.nanmean(phyloP, axis=1)
    for p in positions:
        d[f'phyloP{p}_avg'] = phyloP_mean[p-1]

    # Features: conservation scores histograms for each codon position - phastCons
    phastCons_bins = np.concatenate((np.linspace(0, 0.75, 4), np.linspace(0.8, 1.0, 5)), axis=0)
    for p in positions:
        hist, _ = np.histogram(phastCons[p-1, :], phastCons_bins)
        for i, hist_value in enumerate(hist):
            d[f'phastCons{p}_hist_{phastCons_bins[i]}-{phastCons_bins[i+1]}'] = hist_value

    # Features: conservation scores histograms for each codon position - phyloP
    phyloP_bins = np.concatenate((np.array([-14, -1]), np.linspace(0, 3, 4), np.linspace(3.5, 6, 6)), axis=0)
    for p in positions:
        hist, _ = np.histogram(phyloP[p-1, :], phyloP_bins)
        for i, hist_value in enumerate(hist):
            d[f'phyloP{p}_hist_{phyloP_bins[i]}-{phyloP_bins[i+1]}'] = hist_value

    # Features: histogram of avg in each codon
    phastCons_codons_avg = np.nanmean(phastCons, axis=0)
    hist, _ = np.histogram(phastCons_codons_avg, phastCons_bins)
    for i, hist_value in enumerate(hist):
        d[f'phastCons_codons_hist_{phastCons_bins[i]}-{phastCons_bins[i+1]}'] = hist_value

    phyloP_codons_avg0 = np.nanmean(phyloP, axis=0)
    hist, _ = np.histogram(phyloP_codons_avg0, phyloP_bins)
    for i, hist_value in enumerate(hist):
        d[f'phyloP_codons_hist_{phyloP_bins[i]}-{phyloP_bins[i+1]}'] = hist_value

    return d


def sub_matrix_features(sub_list, weigted_sub_list, sub_name):
    if len(sub_list) == 0:
        sub_avg = weigted_sub_avg = sub_postivies = sub_negatives = sub_ratio = 1
    else:
        # Feature: BLOSUM62 average and frequency weighted-average
        sub_avg = sum(sub_list) / float(len(sub_list))
        weigted_sub_avg = sum(weigted_sub_list) / float(len(weigted_sub_list))

        # Feature: BLOSUM62 count of positives and negatives
        sub_postivies = sum(1 for x in sub_list if x > 0)
        sub_negatives = sum(1 for x in sub_list if x < 0)

        # Feature: BLOSUM62 positives/negatives ratio
        if sub_postivies == 0 or sub_negatives == 0:
            sub_ratio = 0
        else:
            sub_ratio = sub_postivies / float(sub_negatives)

    return {
        f'{sub_name}_avg': sub_avg,
        f'{sub_name}_avg_weighted': weigted_sub_avg,
        f'{sub_name}_positive_num': sub_postivies,
        f'{sub_name}_negative_num': sub_negatives,
        f'{sub_name}_ratio': sub_ratio
    }


def SIFT_features(sift_scores_list, weighted_sift_scores_list):
    if len(sift_scores_list) > 0:
        # Feature: SIFT average
        sift_avg = np.mean(sift_scores_list)

        # Feature: weighted (by frequency) SIFT average
        sift_w_avg = np.mean(weighted_sift_scores_list)

        # Feature: SIFT number of deleterious (score <=0.05)
        sift_deleterious_num = sum(1 for x in sift_scores_list if x <= SIFT_THRESHOLD)

        # Feature: SIFT number of tolerated (score > 0.05)
        sift_tolerated_num = sum(1 for x in sift_scores_list if x > SIFT_THRESHOLD)

        # Feature: deleterious/tolerated ratio
        if sift_tolerated_num == 0 or sift_deleterious_num == 0:
            sift_ratio = 0
        else:
            sift_ratio = sift_deleterious_num / float(sift_tolerated_num)

        # Feature: SIFT "majority-decision" (deleterious/tolerated)
        if sift_deleterious_num > sift_tolerated_num:
            sift_majority = sift_codes.SIFT_DELETERIOUS.value
        elif sift_tolerated_num > sift_deleterious_num:
            sift_majority = sift_codes.SIFT_TOLERATED.value
        else:
            sift_majority = sift_codes.SIFT_TIE.value

    else:
        sift_avg = sift_w_avg = -1
        sift_deleterious_num = 0
        sift_tolerated_num = 0
        sift_ratio = 1
        sift_majority = sift_codes.SIFT_TIE.value

    return {
        'sift_avg': sift_avg,
        'sift_avg_weighted': sift_w_avg,
        'sift_deleterious_num': sift_deleterious_num,
        'sift_tolerated_num': sift_tolerated_num,
        'sift_ratio': sift_ratio,
        'sift_majority': sift_majority
    }


def PolyPhen_features(polyphen_scores_list, polyphen_pred_list, weighted_polyphen_scores_list):
    if len(polyphen_scores_list) > 0:
        # Feature: PolyPhen average
        polyphen_avg = np.mean(polyphen_scores_list)

        # Feature: weighted (by frequency) PolyPhen average
        polyphen_w_avg = np.mean(weighted_polyphen_scores_list)

        # Feature: polyPhen number of benign
        polyphen_benign_num = polyphen_pred_list.count("benign")

        # Feature: polyPhen number of possibly_damaging
        polyphen_possibly_num = polyphen_pred_list.count("possibly_damaging")

        # Feature: polyPhen number of probably_damaging
        polyphen_probably_num = polyphen_pred_list.count("probably_damaging")

        # Feature: polyPhen "majority-decision" (benign/possibly_damaging/probably_damaging/unknown)
        if ((polyphen_benign_num > polyphen_probably_num and polyphen_benign_num > polyphen_possibly_num) or
                (polyphen_benign_num > polyphen_probably_num and polyphen_benign_num == polyphen_possibly_num)):
            polyphen_majority = polyphen_codes.POLYPHEN_BENIGN.value

        elif ((polyphen_probably_num > polyphen_benign_num and polyphen_probably_num > polyphen_possibly_num) or
              (polyphen_probably_num > polyphen_benign_num and polyphen_probably_num == polyphen_possibly_num)):
            polyphen_majority = polyphen_codes.POLYPHEN_PROBABLY.value

        elif polyphen_possibly_num > polyphen_benign_num and polyphen_possibly_num > polyphen_probably_num:
            polyphen_majority = polyphen_codes.POLYPHEN_POSSIBLY.value

        elif polyphen_benign_num == polyphen_probably_num == polyphen_possibly_num:
            polyphen_majority = polyphen_codes.PLOYPHEN_EQUAL.value

        else:
            polyphen_majority = polyphen_codes.POLYPHEN_UNKNOWN.value

    else:
        polyphen_avg = polyphen_w_avg = -1
        polyphen_benign_num = 0
        polyphen_possibly_num = 0
        polyphen_probably_num = 0
        polyphen_majority = polyphen_codes.POLYPHEN_UNKNOWN.value

    return {
        'polyphen_avg': polyphen_avg,
        'polyphen_avg_weighted': polyphen_w_avg,
        'polyphen_benign_num': polyphen_benign_num,
        'polyphen_possibly_num': polyphen_possibly_num,
        'polyphen_probably_num': polyphen_probably_num,
        'polyphen_majority': polyphen_majority
    }


def ClinVar_scores(clinsig_list, clinsig_af):
    valid_scores = []
    valid_scores_weighted = []

    for i in range(len(clinsig_list)):
        sig = clinsig_list[i]
        sig_list = pd.Series(sig.split("&")).unique().tolist()
        # Skipping
        if "not" in sig_list or "" in sig_list:
            continue

        # Determine the alteration clinvar score
        if len(sig_list) == 1 and sig_list[0] == "pathogenic":
            score = clinvar_codes.CLINVAR_PATHOGENIC.value
        elif len(sig_list) == 1 and sig_list[0] == "benign":
            score = clinvar_codes.CLINVAR_BENIGN.value
        elif len(sig_list) == 2 and "benign" in sig_list and "likely" in sig_list:
            score = clinvar_codes.CLINVAR_LIKELY_BENIGN.value
        elif len(sig_list) == 2 and "pathogenic" in sig_list and "uncertain" in sig_list:
            score = clinvar_codes.CLINVAR_LIKELY_PATHOGENIC.value
        elif len(sig_list) == 2 and "pathogenic" in sig_list and "other" in sig_list:
            score = clinvar_codes.CLINVAR_PATHOGENIC_OTHER.value
        else:
            score = clinvar_codes.CLINVAR_UNCERTAIN.value  # value of 0

        valid_scores.append(score)
        score_af = clinsig_af[i]
        valid_scores_weighted.append(score * score_af)

    # ===Feature: Avg. and weighted avg. ClinVar score===#
    if len(valid_scores) == 0:
        avg_clinvar_score = 0
        avg_w_clinvar_score = 0
    else:
        avg_clinvar_score = np.mean(valid_scores)
        avg_w_clinvar_score = np.mean(valid_scores_weighted)

    return {
        'avg_clinvar_score': avg_clinvar_score,
        'avg_clinvar_weighted': avg_w_clinvar_score
    }


def entropy_features(maf_list):
    # Calculates a normalized Shannon entropy (from Miller et al, 2015) of nonsyn SNPs distributed across instances
    maf = np.array(maf_list)
    if np.sum(maf) == 0:
        # if no SNPs- each instance has the prob. = max. entropy ln(n)
        e = np.log(len(maf))
    else:
        # Filter out nans; scipy.stats.entropy automatically normalizes the input
        # We divide the result by ln(|x|) to account for different sized inputs
        e = entropy(maf[~np.isnan(maf)]) / np.log(len(maf))

    return {'snp_nonsyn_entropy': e}


def pseudo_dNdS_features(ref_seq, Nd, Sd):
    N, S = seq_ns(ref_seq)  # Reference expected syn/nonsyn per site
    PN = 0 if N == 0 else Nd / float(N)  # Proportion of nonsyn
    PS = 0 if S == 0 else Sd / float(S)  # Proportion of syn

    # num of nonsyn substitutions per nonsyn site
    dN = -0.75 * (np.log(1 - 4 * PN / float(3)))

    # num of syn substitutions per syn site
    if 4 * PS / float(3) >= 1:
        dS = 1
    else:
        dS = -0.75 * (np.log(1 - 4 * PS / float(3)))

    if dN == 0 or dS == 0:
        dN_dS = 1  # There isn't enough information to calculate dN/dS (1 is a neutral value)
    else:
        dN_dS = dN / dS
        if dN_dS == np.nan:
            dN_dS = 1  # There isn't enough information to calculate dN/dS (1 is a neutral value)

    return {
        'pseudo_nonsyn': dN,
        'pseudo_syn': dS,
        'pseudo_dNdS': dN_dS
    }


def pfam_emission_prob_features(hmm_prob_dict, state):
    # Max. emission probability + emission prob. for each amino acid
    probs = hmm_prob_dict[state]
    d = {'pfam_prob_max': max(probs)}
    d.update({f'pfam_prob_{pfam_aa_order[i]}': prob for i, prob in enumerate(probs)})
    return d


def pfam_conserved_state_feature(state, con_states_dict):
    # is state is conserved according to Pfam?
    return {'is_pfam_conserved': state in con_states_dict}


def instance_individuals_100way_change_features(maf_list, aa_ref_hist, jsd100way_list):

    # Computing Orthologus conservation in different ways (from 100way-ucsc alignment)
    # Computing Paralogus conservartion in different ways (from different instances)
    # Combining both to measurments that maximize ortho. con. and minimize para. con.

    ## Paralogus ##

    # fraction of change across instances

    # determine majority aa (index of one of the majority)
    minor_counts = 0
    max_pos = aa_ref_hist.index(max(aa_ref_hist))
    for i in range(len(aa_ref_hist)):
        if i == max_pos:
            continue
        minor_counts += aa_ref_hist[i]

    instances_change_frac = minor_counts / float(np.sum(aa_ref_hist))

    # Feature: entropy of ref AA
    aa_ref_entropy = SE_hist(aa_ref_hist)

    # JSD of ref AA
    aa_ref_jsd = JSD_hist(aa_ref_hist, background=JSD_background.BLOSUM62)

    ## Orthologus ##

    # first remove -1 illegal scores of JSD mismatch (positions where JSD alignment didn't match, I added -1):
    jsd100way_list_no_mismatch = [i for i in jsd100way_list if i != -1]

    # median JSD score across 100way vertbrates
    med_jsd = 0 if len(jsd100way_list_no_mismatch) == 0 else np.median(jsd100way_list_no_mismatch)

    # Histogram of JSD score across 100way vertebrates
    jsd_median_bins = [0, 0.5, 0.6, 0.7, 0.8, 1]
    jsd_median_hist = np.histogram(jsd100way_list_no_mismatch, bins=jsd_median_bins)[0]

    # Functional measurements of both
    # ratio: change across instances / change across individuals(MAF)

    # low MAF (orthologues), high instances change (paralogous) = SDPs
    if np.sum(maf_list) == 0:
        avg_maf_overall = 0.0000001  # set the minimal non-zero in our data
    else:
        avg_maf_overall = np.sum(maf_list) / float(len(maf_list))

    max_entropy = SE_hist([0] * len(AMINO_ACIDS))

    return dict(
        [
            ('instances_change_frac', instances_change_frac),
            ('aa_ref_SE', aa_ref_entropy),
            ('aa_ref_jsd', aa_ref_jsd),
            ('med_jsd_100way_blosum', med_jsd)
        ] +
        [
            (f'jsd_median_hist_{jsd_median_bins[i]}-{jsd_median_bins[i + 1]}', jsd_median_hist[i])
            for i in range(len(jsd_median_bins) - 1)
        ] +
        [
            ('instances_individuals_change_ratio', instances_change_frac / float(avg_maf_overall)),

            # high JSD (orthologues), high instances change (paralogous) = SDPs
            # we want high MAF -> small 1 - MAF, high JSD
            ('jsd_100way_instances_major_ratio', med_jsd / float(1 - instances_change_frac)),

            # high JSD (orthologues), high shannon entropy (paralogous) = SDPs
            ('jsd_mul_aa_ref_SE', med_jsd * aa_ref_entropy),

            # high JSD (orthologues), low diff. of max SE to shannon entropy (paralogous) = SDPs
            ('jsd_SE_diff_ratio', med_jsd / float(max_entropy - aa_ref_entropy)),

            # high JSD (orthologues), high shannon entropy (paralogous) = SDPs
            ('jsd_SE_sum', med_jsd + (aa_ref_entropy / float(max_entropy))),

            # high shannon entropy (paralogous), low diff. of max JSD to avg JSD (orthologues) = SDPs
            ('SE_jsd_diff_ratio', aa_ref_entropy / float(1 - med_jsd)),

            # high JSD (orthologues), low JSD (paralogoues) = SDPs
            ('jsds_ratio', med_jsd / float(aa_ref_jsd)),

            # high difference between orthoulogus (more conserved) and paralogous (less conserved)
            ('jsds_subtraction', med_jsd - aa_ref_jsd)
        ]

    )


def aa_identity_features(aa_ref_hist, type_str):
    # aa identity histogram and probability
    if np.sum(aa_ref_hist) == 0:
        aa_ref_prob = aa_ref_hist
    else:
        aa_ref_prob = np.asarray(aa_ref_hist) / float(np.sum(aa_ref_hist))

    d = {f'{type_str}_hist_{aa}': aa_ref_hist_i for aa, aa_ref_hist_i in zip(AMINO_ACIDS, aa_ref_hist)}
    d.update({f'{type_str}_prob_{aa}': aa_ref_prob_i for aa, aa_ref_prob_i in zip(AMINO_ACIDS, aa_ref_prob)})
    return d


def major_allele_charge(aa_ref_hist):
    # ===Feature: major allele aa charge counts===#
    charge_positive_count = charge_negative_count = charge_neutral_count = 0
    for i in range(len(AMINO_ACIDS)):
        aa_count = aa_ref_hist[i]
        if aa_count > 0:
            charge = aa_charge_dict[AMINO_ACIDS[i]]
            if charge.value == 0:
                charge_neutral_count += aa_count
            elif charge.value == 1:
                charge_positive_count += aa_count
            else:
                charge_negative_count += aa_count

    # ===Feature: major allele majority charge===#
    charge_majority = aa_charge.NEUTRAL.value
    if charge_positive_count > charge_neutral_count and charge_positive_count > charge_negative_count:
        charge_majority = aa_charge.POSITIVE.value
    elif charge_negative_count > charge_neutral_count and charge_negative_count > charge_positive_count:
        charge_majority = aa_charge.NEGATIVE.value

    return {
        'aa_ref_charge_positive_count': charge_positive_count,
        'aa_ref_charge_negative_count': charge_negative_count,
        'aa_ref_charge_neutral_count': charge_neutral_count,
        'aa_ref_charge_majority': charge_majority
    }


def major_allele_functional_group(aa_ref_hist):
    # major allele aa functional group counts
    func_counters = [0] * (len(aa_functional_group) - 1)  # Major allele is never a stop codon
    for i in range(len(AMINO_ACIDS)):
        aa_count = aa_ref_hist[i]
        if aa_count > 0:
            func_group_num = aa_functional_group_dict[
                AMINO_ACIDS[i]].value  # getting numeric functional group value
            if func_group_num == aa_functional_group.STOP.value:  # Major allele is never a stop codon
                continue
            func_counters[func_group_num] += aa_count

    return {
        k: v for k, v in zip(
            [f'aa_ref_{group}_count' for group in aa_functional_group if group != aa_functional_group.STOP],
            func_counters
        )
    }


def sub_diff_functional_group(ref_alt_pairs):
    # ===Features: count and frequency staying in functional group Vs. moving to other group===#
    stay_cnt = stay_cnt_freq = move_cnt = move_cnt_freq = 0

    for (ref, alt, af) in ref_alt_pairs:
        ref_func_group = aa_functional_group_dict[ref].value
        alt_func_group = aa_functional_group_dict[alt].value
        if ref_func_group == alt_func_group:
            stay_cnt += 1
            stay_cnt_freq += af
        else:
            move_cnt += 1
            move_cnt_freq += af

    # ===Features: functional groups transitions counts===#
    transitions_vec_size = (len(aa_functional_group) - 1) * len(
        aa_functional_group)  # excluding transitions from STOP codons
    transitions_vec = [0] * transitions_vec_size

    for (ref, alt, af) in ref_alt_pairs:
        ref_func_group = aa_functional_group_dict[ref].value
        alt_func_group = aa_functional_group_dict[alt].value
        # Calculate counter position on the vector (ref_func_group is never STOP = 5)
        trans_vec_i = ref_func_group * (len(aa_functional_group) - 1)
        trans_vec_i += alt_func_group
        transitions_vec[trans_vec_i] += 1

    keys = [f'sub_func_group_trans_{i}-{j}' for i in range(len(aa_functional_group) - 1) for j in range(len(aa_functional_group))]

    return dict([
        ('sub_func_group_stay_cnt', stay_cnt),
        ('sub_func_group_stay_freq', stay_cnt_freq),
        ('sub_func_group_move_cnt', move_cnt),
        ('sub_func_group_move_freq', move_cnt_freq)
    ] + [(k, v) for k, v in zip(keys, transitions_vec)])


def major_allele_hydrophobicity(aa_ref_hist):
    # major allele hydrophicity average, hydrophobic and polar counts
    h_sum = h_cnt = hydrophobic_cnt = polar_charge_cnt = 0
    for i in range(len(AMINO_ACIDS)):
        aa_count = aa_ref_hist[i]
        if aa_count > 0:
            hindex = kd.get(AMINO_ACIDS[i], 0)
            h_sum += hindex * aa_count
            h_cnt += aa_count

            if hindex > 0:
                hydrophobic_cnt += aa_count
            else:
                polar_charge_cnt += aa_count

    h_avg = 0 if h_cnt == 0 else h_sum / float(h_cnt)

    return {
        'hindex_avg': h_avg,
        'hindex_pos_cnt': hydrophobic_cnt,
        'hindex_neg_cnt': polar_charge_cnt
    }


def sub_diff_hydrophobicity(ref_alt_pairs):
    # hydrophicity difference average and weighted average
    hindex_diff_sum = hindex_diff_sum_weighted = hindex_diff_cnt = 0
    for (ref, alt, af) in ref_alt_pairs:
        hindex_diff = kd.get(alt, 0) - kd.get(ref, 0)
        hindex_diff_sum += hindex_diff
        hindex_diff_sum_weighted += hindex_diff * af
        hindex_diff_cnt += 1

    if hindex_diff_cnt == 0:
        hindex_diff_avg = hindex_diff_avg_weighted = 0
    else:
        hindex_diff_avg = hindex_diff_sum / float(hindex_diff_cnt)
        hindex_diff_avg_weighted = hindex_diff_sum_weighted / float(hindex_diff_cnt)

    return {
        'sub_diff_hindex_avg': hindex_diff_avg,
        'sub_diff_hindex_avg_weighted': hindex_diff_avg_weighted
    }


def major_allele_volume(aa_ref_hist):
    # major allele volume average, tiny, small and big counts
    vol_sum = vol_cnt = tiny_cnt = small_cnt = big_cnt = 0
    for i in range(len(AMINO_ACIDS)):
        aa_count = aa_ref_hist[i]
        if aa_count > 0:
            volume = aa_volume[AMINO_ACIDS[i]]
            vol_sum += volume * aa_count
            vol_cnt += aa_count

            vol_group = aa_volume_group_dict[AMINO_ACIDS[i]]
            if vol_group == aa_volume_group.TINY:
                tiny_cnt += aa_count
            elif vol_group == aa_volume_group.SMALL:
                small_cnt += aa_count
            elif vol_group == aa_volume_group.BIG:
                big_cnt += aa_count

    vol_avg = 0 if vol_cnt == 0 else vol_sum / float(vol_cnt)

    return {
        'vol_avg': vol_avg,
        'vol_tiny_cnt': tiny_cnt,
        'vol_small_cnt': small_cnt,
        'vol_big_cnt': big_cnt
    }


def sub_diff_volume(ref_alt_pairs):
    # ===Feature: volume difference average and weighted average===#
    volume_diff_sum = 0
    volume_diff_sum_weighted = 0
    volume_diff_cnt = 0
    for (ref, alt, af) in ref_alt_pairs:
        ref_vol = aa_volume[ref]
        alt_vol = aa_volume[alt]
        vol_diff = (ref_vol - alt_vol)
        volume_diff_sum += vol_diff
        volume_diff_sum_weighted += vol_diff * af
        volume_diff_cnt += 1

    if volume_diff_cnt == 0:
        volume_diff_avg = volume_diff_avg_weighted = 0
    else:
        volume_diff_avg = volume_diff_sum / float(volume_diff_cnt)
        volume_diff_avg_weighted = volume_diff_sum_weighted / float(volume_diff_cnt)

    return {
        'sub_diff_vol_avg': volume_diff_avg,
        'sub_diff_vol_avg_weighted': volume_diff_avg_weighted
    }


def major_allele_propensity(aa_ref_hist):
    prop_sum = [0, 0, 0]
    prop_cnt = 0
    prop_majority_counts = [0, 0, 0]
    for i in range(len(AMINO_ACIDS)):
        aa_count = aa_ref_hist[i]
        if aa_count > 0:
            curr_prop = propensity_chou_fasman[AMINO_ACIDS[i]]
            mul_curr_prop = [x * aa_count for x in curr_prop]
            prop_sum = [sum(x) for x in zip(prop_sum, mul_curr_prop)]
            prop_cnt += aa_count

            if curr_prop[aa_propensity.ALPHA_HELIX.value] == max(curr_prop):
                prop_majority_counts[aa_propensity.ALPHA_HELIX.value] += 1
            if curr_prop[aa_propensity.BETA_SHEET.value] == max(curr_prop):
                prop_majority_counts[aa_propensity.BETA_SHEET.value] += 1
            if curr_prop[aa_propensity.TURN.value] == max(curr_prop):
                prop_majority_counts[aa_propensity.TURN.value] += 1

    # ===Feature: major allele propensity avgs===#
    if prop_cnt == 0:
        prop_avg = [0, 0, 0]
    else:
        prop_avg = [x / float(prop_cnt) for x in prop_sum]

    # ===Feature: major allele majority propensity===#
    max_idx = np.where(np.array(prop_majority_counts) == max(prop_majority_counts))[0]
    majority_vec = [0, 0, 0]
    for i in max_idx:
        majority_vec[i] = 1  # put 1 in the propensities that has max. count

    return {
        'aa_ref_alpha_prop_avg': prop_avg[0],
        'aa_ref_beta_prop_avg': prop_avg[1],
        'aa_ref_turn_prop_avg': prop_avg[2],
        'aa_ref_alpha_is_majority': majority_vec[0],
        'aa_ref_beta_is_majority': majority_vec[1],
        'aa_ref_turn_is_majority': majority_vec[2]
    }


def sub_diff_propensity(ref_alt_pairs):
    # ===Feature: propensity difference average===#
    prop_vec_sum = [0, 0, 0]
    prop_vec_sum_weighted = [0, 0, 0]
    prop_cnt = 0
    for (ref, alt, af) in ref_alt_pairs:
        ref_struct = propensity_chou_fasman[ref]
        alt_struct = propensity_chou_fasman[alt]
        prop_diff = [(x - y) for (x, y) in zip(ref_struct, alt_struct)]
        prop_diff_weighted = [(x - y) * af for (x, y) in zip(ref_struct, alt_struct)]
        prop_vec_sum = [(x + y) for (x, y) in zip(prop_vec_sum, prop_diff)]
        prop_vec_sum_weighted = [(x + y) for (x, y) in zip(prop_vec_sum_weighted, prop_diff_weighted)]

        prop_cnt += 1

    if prop_cnt == 0:
        prop_vec_avg = prop_vec_avg_weighted = [0, 0, 0]
    else:
        prop_vec_avg = [(x / float(prop_cnt)) for x in prop_vec_sum]
        prop_vec_avg_weighted = [(x / float(prop_cnt)) for x in prop_vec_sum_weighted]

    return {
        'sub_diff_prop_avg_alpha': prop_vec_avg[0],
        'sub_diff_prop_avg_beta': prop_vec_avg[1],
        'sub_diff_prop_avg_turn': prop_vec_avg[2],
        'sub_diff_prop_avg_alpha_weighed': prop_vec_avg_weighted[0],
        'sub_diff_prop_avg_beta_weighed': prop_vec_avg_weighted[1],
        'sub_diff_prop_avg_turn_weighed': prop_vec_avg_weighted[2]
    }


def major_allele_h_bonds(aa_ref_hist):
    # avg donor and acceptor H-bond potential
    donor_sum = acceptor_sum = bonds_cnt = 0
    for i in range(len(AMINO_ACIDS)):
        aa_count = aa_ref_hist[i]
        if aa_count > 0:
            donor_sum += (aa_h_bond_donor[AMINO_ACIDS[i]] * aa_count)
            acceptor_sum += (aa_h_bond_acceptor[AMINO_ACIDS[i]] * aa_count)
            bonds_cnt += aa_count

    if bonds_cnt == 0:
        donor_avg = 0
        acceptor_avg = 0
    else:
        donor_avg = donor_sum / float(bonds_cnt)
        acceptor_avg = acceptor_sum / float(bonds_cnt)

    return {
        'H_bond_donor_avg': donor_avg,
        'H_bond_acceptor_avg': acceptor_avg
    }


def sub_diff_h_bonds(ref_alt_pairs):
    # acceptor and donor diff average and weighted average
    donor_diff_sum = donor_diff_sum_weighted = acceptor_diff_sum = acceptor_diff_sum_weighted = diff_cnt = 0
    for (ref, alt, af) in ref_alt_pairs:
        donor_diff = aa_h_bond_donor[ref] - aa_h_bond_donor[alt]
        donor_diff_sum += donor_diff
        donor_diff_sum_weighted += donor_diff * af

        ref_acceptor = aa_h_bond_acceptor[ref]
        alt_acceptor = aa_h_bond_acceptor[alt]
        acceptor_diff = (ref_acceptor - alt_acceptor)
        acceptor_diff_sum += acceptor_diff
        acceptor_diff_sum += acceptor_diff * af

        diff_cnt += 1

    if diff_cnt == 0:
        donor_diff_avg = donor_diff_avg_weighted = 0
        acceptor_diff_avg = acceptor_diff_avg_weighted = 0
    else:
        donor_diff_avg = donor_diff_sum / float(diff_cnt)
        donor_diff_avg_weighted = donor_diff_sum_weighted / float(diff_cnt)
        acceptor_diff_avg = acceptor_diff_sum / float(diff_cnt)
        acceptor_diff_avg_weighted = acceptor_diff_sum_weighted / float(diff_cnt)

    return {
        'donor_diff_avg': donor_diff_avg,
        'donor_diff_avg_weighted': donor_diff_avg_weighted,
        'acceptor_diff_avg': acceptor_diff_avg,
        'acceptor_diff_avg_weighted': acceptor_diff_avg_weighted
    }


def spider_solvent_acc_pred(spider_dict):
    # Accessible Surface Area (solvent accessibility) mean/std
    return {
        'solvent_acc_avg': np.nanmean(spider_dict["spider2-ASA"]),
        'solvent_acc_std': np.nanstd(spider_dict["spider2-ASA"])
    }


def spider_contact_number_pred(spider_dict):
    return {
        'hsa2_cn_avg': np.nanmean(spider_dict["spider2-hsa2_CN"]),  # contact number for Cα-Cα mean
        'hsa2_cn_std': np.nanstd(spider_dict["spider2-hsa2_CN"]),   # contact number for Cα-Cα std
        'hsb2_cn_avg': np.nanmean(spider_dict["spider2-hsb2_CN"]),  # contact number for Cα-Cβ mean
        'hsb2_cn_std': np.nanstd(spider_dict["spider2-hsb2_CN"])    # contact number for Cα-Cβ std
    }


def spider_angles_pred(spider_dict):
    return {
        'backbone_Phi_angle_avg': np.nanmean(spider_dict["spider2-angle_Phi"]),     # backbone Phi angle mean
        'backbone_Phi_angle_std': np.nanstd(spider_dict["spider2-angle_Phi"]),      # backbone Phi angle std
        'backbone_Psi_angle_avg': np.nanmean(spider_dict["spider2-angle_Psi"]),     # backbone Psi angle mean
        'backbone_Psi_angle_std': np.nanstd(spider_dict["spider2-angle_Psi"]),      # backbone Psi angle std
        'c-alpha_tau_angle_avg': np.nanmean(spider_dict["spider2-angle_tau"]),      # c-alpha angle (i-2=>i+1) mean
        'c-alph_tau_angle_std': np.nanstd(spider_dict["spider2-angle_tau"]),        # c-alpha angle (i-2=>i+1) std
        'c-alpha_theta_angle_avg': np.nanmean(spider_dict["spider2-angle_theta"]),  # c-alpha angle (i-1=>i+1) mean
        'c-alph_theta_angle_std': np.nanstd(spider_dict["spider2-angle_theta"])     # c-alpha angle (i-1=>i+1) std
    }


def spider_struct_pred(spider_dict):

    # major allele majority propensity
    values, counts = np.unique(np.array(spider_dict["spider2-2nd_struct"]), return_counts=True)

    # Note: An earlier iteration of the code returned all 3 spd_*_is_majority keys as 1s
    # in the absence of any spider 2nd struct infomation, while at the same time returning all the
    # *_prob_avg/std keys as nans
    # We follow the same logic here by setting the major allele (which would normally be one of H/E/C as HEC
    # in case of missing information.
    # Note also that we cannot use np.argmax since it only returns a 'single' index
    maj_allele = 'HEC' if len(counts) == 0 else values[np.where(counts == np.max(counts))]

    return {
        'helix_prob_avg': np.nanmean(spider_dict["spider2-helix_prob"]),  # helix prob. mean
        'helix_prob_std': np.nanstd(spider_dict["spider2-helix_prob"]),   # helix prob. std
        'sheet_prob_avg': np.nanmean(spider_dict["spider2-sheet_prob"]),  # sheet prob. mean
        'sheet_prob_std': np.nanstd(spider_dict["spider2-sheet_prob"]),   # sheet prob. std
        'turn_prob_avg': np.nanmean(spider_dict["spider2-turn_prob"]),    # turn prob. mean
        'turn_prob_std': np.nanstd(spider_dict["spider2-turn_prob"]),     # turn prob. std
        'spd_helix_is_majority': int('H' in maj_allele),                  # 0/1
        'spd_sheet_is_majority': int('E' in maj_allele),                  # 0/1
        'spd_turn_is_majority': int('C' in maj_allele)                    # 0/1
    }


def spider_half_sphere_exposure_pred(spider_dict):
    return {
        'hsa2_HSE-up_avg': np.mean(spider_dict["spider2-hsa2_HSEu"]),
        'hsa2_HSE-up_std': np.std(spider_dict["spider2-hsa2_HSEu"]),
        'hsa2_HSE-down_avg': np.mean(spider_dict["spider2-hsa2_HSEu"]),  # TODO: should be HSEd, bug?
        'hsa2_HSE-down_std': np.std(spider_dict["spider2-hsa2_HSEd"]),
        'hsb2_HSE-up_avg': np.mean(spider_dict["spider2-hsb2_HSEu"]),
        'hsb2_HSE-up_std': np.std(spider_dict["spider2-hsb2_HSEu"]),
        'hsb2_HSE-down_avg': np.mean(spider_dict["spider2-hsb2_HSEd"]),
        'hsb2_HSE-down_std': np.std(spider_dict["spider2-hsb2_HSEd"])
    }


def whole_domain_conservation(states_dict):
    # phastCons and PhyloP whole-domain mean/std
    return {
        'whole_domain_phastCons_avg': states_dict["_phastCons_mean"],
        'whole_domain_phastCons_std': states_dict["_phastCons_std"],
        'whole_domain_phyloP_avg': states_dict["_phyloP_mean"],
        'whole_domain_phyloP_std': states_dict["_phyloP_std"]
    }


def domain_location_features(state, max_state):

    # ===Feature: the location in the domain: beginning/middle/end===#
    location_list = [0, 0, 0]
    BEGIN_POS = 0
    MIDDLE_POS = 1
    END_POS = 2
    domain_location_bins = np.histogram(np.arange(1, max_state), bins=3)[1]
    if state < domain_location_bins[1]:
        location_list[BEGIN_POS] = 1
    elif state > domain_location_bins[2]:
        location_list[END_POS] = 1
    else:
        location_list[MIDDLE_POS] = 1

    return {
        'domain_pos': state,
        'domain_length': max_state,
        'domain_pos_location_begin': location_list[0],
        'domain_pos_location_middle': location_list[1],
        'domain_pos_location_end': location_list[2]
    }


def protein_location_features(protein_pos_list, protein_len_list):
    # Avg. protein total length, counts of the location in the protein: beginning/middle/end
    location_list = [0, 0, 0]
    BEGIN_POS = 0
    MIDDLE_POS = 1
    END_POS = 2
    for i in range(len(protein_pos_list)):
        prot_location_bins = np.histogram(np.arange(1, protein_len_list[i]), bins=3)[1]
        if protein_pos_list[i] < prot_location_bins[1]:
            location_list[BEGIN_POS] += 1
        elif protein_pos_list[i] > prot_location_bins[2]:
            location_list[END_POS] += 1
        else:
            location_list[MIDDLE_POS] += 1

    # Normalize to ratios
    location_list_norm = np.array(location_list) / sum(location_list)

    return {
        'prot_avg_length': np.mean(protein_len_list),
        'prot_pos_location_begin': location_list_norm[0],
        'prot_pos_location_middle': location_list_norm[1],
        'prot_pos_location_end': location_list_norm[2]
    }


if __name__ == '__main__':

    with open(os.path.join(os.path.dirname(data.__file__), 'BLOSUM62_dict.pik'), 'rb') as f:
        blosum62_dict = pickle.load(f)

    with open(os.path.join(os.path.dirname(data.__file__), 'PAM40_dict.pik'), 'rb') as f:
        pam40_dict = pickle.load(f)

    with open(PROB_DICT, 'rb') as f:
        prob_dict = pickle.load(f)

    features_list = []
    an_str = POPULATIONS_ANS
    ac_str = POPULATIONS_ACS

    for pik in glob.glob(f'{HMM_STATES_FOLDER}/*.pik'):
        domain = os.path.splitext(os.path.basename(pik))[0]


        hmm_prob_dict = prob_dict[domain]

        with open(pik, 'rb') as handle:
            states_dict = pickle.load(handle)

        # Create af_adj flat dict
        states_af_adj_dict = defaultdict(list)
        for state in states_dict.keys():
            if str(state).startswith('_'): continue
            for d in states_dict[state]:
                states_af_adj_dict[state].append(d["af_adj"])

        # scale the af_dict
        states_MAF_adj_dict_scaled = defaultdict(list)
        for state in states_dict.keys():
            if str(state).startswith('_'): continue
            state_len = len(states_dict[state])
            for d in states_dict[state]:
                states_MAF_adj_dict_scaled[state].append(float(d["af_adj"] / state_len))

        # Create a dict of conserved states
        con_states_dict = {}
        con_threshold = 0.5
        for state in hmm_prob_dict:
            prob_list = hmm_prob_dict[state]
            for i in range(len(prob_list)):
                p = prob_list[i]
                if p > con_threshold:
                    major_allele = pfam_aa_order[i]
                    con_states_dict[state] = major_allele

        # Adding states features
        for state in states_dict:

            features_dict = {}

            if str(state).startswith('_'): continue
            state_id = domain + "_" + str(state)

            # Init counters & paramters
            maf_list = []
            sites_aa_alter_num = 0
            sites_snp_alter_num = 0
            sites_aa_num = len(states_dict[state])
            sites_snp_num = 3 * sites_aa_num
            sites_poly_aa_num = 0  # The number of different aa in all the altered sites (most are 1)
            sites_poly_aa_several = 0

            # Rare-poly-counters
            rare_5_num = 0
            rare_05_num = 0
            rare_005_num = 0

            # Conservation params
            phastCons_dict = defaultdict(list)
            phyloP_dict = defaultdict(list)
            jsd100way_list = []

            # SPIDER params
            spider_dict = defaultdict(list)

            # BLOSUM62_params
            blosum62_list = []
            weigted_blosum62_list = []

            # PAM40_params
            pam40_list = []
            weigted_pam40_list = []

            # dn/ds counters and variables
            ref_seq = ""
            Nd = 0
            Sd = 0

            # SIFT params
            sift_scores_list = []
            weighted_sift_scores_list = []

            # PolyPhen params
            polyphen_scores_list = []
            weighted_polyphen_scores_list = []
            polyphen_pred_list = []

            # clinVar params
            clinsig_list = []
            clinsig_af = []

            # Major allele params
            aa_ref_hist = [0] * len(AMINO_ACIDS)

            # Substitution params
            aa_alt_hist = [0] * len(AMINO_ACIDS)
            aa_alt_prob = [0] * len(AMINO_ACIDS)
            aa_alt_prob_avg = [0] * len(AMINO_ACIDS)
            ref_alt_pairs = []

            # protein position params
            protein_pos_list = []
            protein_len_list = []

            # Populations variables
            ac_sum = [0] * len(ac_str)
            ac_sum_syn = [0] * len(ac_str)
            ac_sum_nonsyn = [0] * len(ac_str)
            an_list = [[] for i in range(len(an_str))]
            pop_maf_list = [[] for i in range(len(an_str))]
            pop_maf_syn_list = [[] for i in range(len(an_str))]
            pop_maf_nonsyn_list = [[] for i in range(len(an_str))]

            # Iterating the state dict to get properties
            for d in states_dict[state]:

                # a list of all maf per instance
                maf_list.append(d["af_adj"])

                # Creating a position pseudo-ref sequence
                ref_codon = d["bp_ref"]
                ref_seq = ref_seq + ref_codon

                # Calculating frequency-based N/S
                bp_af_adj_dict = d["bp_af_adj_dict"]
                for alt_codon in bp_af_adj_dict.keys():
                    alt_aa = codon_table[alt_codon]
                    # syn
                    if alt_aa == d["aa_ref"]:
                        Sd += bp_af_adj_dict[alt_codon]
                    # Non-syn
                    else:
                        Nd += bp_af_adj_dict[alt_codon]

                # Major allele parameters
                aa_ref = d["aa_ref"]
                aa_ref_pos = AMINO_ACIDS.index(aa_ref)
                aa_ref_hist[aa_ref_pos] += 1

                # Conservation scores
                phastCons_curr_list = d["phastCons"]
                if len(phastCons_curr_list) > 0:
                    phastCons_dict[1].append(phastCons_curr_list[0])
                if len(phastCons_curr_list) > 1:
                    phastCons_dict[2].append(phastCons_curr_list[1])
                else:
                    phastCons_dict[2].append(np.nan)
                if len(phastCons_curr_list) > 2:
                    phastCons_dict[3].append(phastCons_curr_list[2])
                else:
                    phastCons_dict[3].append(np.nan)

                phyloP_curr_list = d["phyloP"]
                if len(phyloP_curr_list) > 0:
                    phyloP_dict[1].append(phyloP_curr_list[0])
                if len(phyloP_curr_list) > 1:
                    phyloP_dict[2].append(phyloP_curr_list[1])
                else:
                    phyloP_dict[2].append(np.nan)
                if len(phyloP_curr_list) > 2:
                    phyloP_dict[3].append(phyloP_curr_list[2])
                else:
                    phyloP_dict[3].append(np.nan)

                jsd100way_list.append(d["100-way-BLOSUM_JSD"])

                # SPIDER parameters (add only if exist)
                if "spider2-2nd_struct" in d:
                    spider_dict["spider2-2nd_struct"].append(d["spider2-2nd_struct"])
                    spider_dict["spider2-helix_prob"].append(float(d["spider2-helix_prob"]))
                    spider_dict["spider2-sheet_prob"].append(float(d["spider2-sheet_prob"]))
                    spider_dict["spider2-turn_prob"].append(float(d["spider2-turn_prob"]))
                    spider_dict["spider2-angle_Phi"].append(float(d["spider2-angle_Phi"]))
                    spider_dict["spider2-angle_Psi"].append(float(d["spider2-angle_Psi"]))
                    spider_dict["spider2-angle_tau"].append(float(d["spider2-angle_tau"]))
                    spider_dict["spider2-angle_theta"].append(float(d["spider2-angle_theta"]))
                    spider_dict["spider2-ASA"].append(float(d["spider2-ASA"]))
                    spider_dict["spider2-hsa2_HSEu"].append(float(d["spider2-hsa2_HSEu"]))
                    spider_dict["spider2-hsa2_HSEd"].append(float(d["spider2-hsa2_HSEd"]))
                    spider_dict["spider2-hsb2_HSEu"].append(float(d["spider2-hsb2_HSEu"]))
                    spider_dict["spider2-hsb2_HSEd"].append(float(d["spider2-hsb2_HSEd"]))
                    spider_dict["spider2-hsa2_CN"].append(float(d["spider2-hsa2_CN"]))
                    spider_dict["spider2-hsb2_CN"].append(float(d["spider2-hsb2_CN"]))

                protein_pos_list.append(d["prot_pos"])
                protein_len_list.append(d["prot_len"])

                if d["af_adj"] > 0:
                    sites_aa_alter_num += 1
                    sites_snp_alter_num += len(d["an_adj"])

                    # Number of different polymorphisms at this site
                    site_poly_num = len(d["alterations_af_adj_dict"].keys())
                    sites_poly_aa_num += site_poly_num
                    if site_poly_num > 1:
                        sites_poly_aa_several += 1

                    # Rare poly features

                    for alt_codon in bp_af_adj_dict.keys():
                        # Add to counters only nonsyn SNPs
                        if codon_table[alt_codon] != codon_table[ref_codon]:
                            if bp_af_adj_dict[alt_codon] < MAFT_005:
                                rare_005_num += 1
                                rare_05_num += 1
                                rare_5_num += 1
                            elif bp_af_adj_dict[alt_codon] < MAFT_05:
                                rare_05_num += 1
                                rare_5_num += 1
                            elif bp_af_adj_dict[alt_codon] < MAFT_5:
                                rare_5_num += 1

                    # Alt, BLOSUM62 and PAM40 features
                    ref = d["aa_ref"]
                    for alt in d["alterations_af_adj_dict"].keys():
                        af_adj = np.mean(d["alterations_af_adj_dict"][alt])
                        # BLOSUM
                        blosum_val = blosum62_dict[ref][alt]
                        blosum62_list.append(blosum_val)
                        weigted_blosum62_list.append(blosum_val * af_adj)
                        # PAM
                        pam_val = pam40_dict[ref][alt]
                        pam40_list.append(pam_val)
                        weigted_pam40_list.append(pam_val * af_adj)
                        # Alt aa counts
                        aa_alt_pos = AMINO_ACIDS.index(alt)
                        aa_alt_hist[aa_alt_pos] += 1
                        # Alt aa prob.
                        aa_alt_prob[aa_alt_pos] += af_adj
                        # ref-alt pairs
                        ref_alt_pairs.append((ref, alt, af_adj))

                    # SIFT
                    sift_list = d["SIFT"]
                    for i in range(len(sift_list)):
                        s = sift_list[i]
                        if s != "":
                            try:
                                s_af = bp_af_adj_dict[d["bp_list"][i]]
                            except:
                                # The major allele was replaced, no score available for the correct substitution
                                continue
                            sift_score = float(s[s.find("(") + 1:s.find(")")])
                            sift_scores_list.append(sift_score)
                            weighted_sift_scores_list.append(sift_score * s_af)

                    # PolyPhen
                    polyphen_list = d["PolyPhen"]
                    for i in range(len(polyphen_list)):
                        s = polyphen_list[i]
                        if s != "":
                            try:
                                s_af = bp_af_adj_dict[d["bp_list"][i]]
                            except:
                                # The major allele was replaced, no score available for the correct substitution
                                continue
                            polyphen_score = float(s[s.find("(") + 1:s.find(")")])
                            polyphen_scores_list.append(polyphen_score)
                            weighted_polyphen_scores_list.append(polyphen_score * s_af)
                            polyphen_pred_list.append(s[:s.find("(")])

                    # clinVar
                    curr_clinsig_list = d["clin_sig"]
                    for i in range(len(curr_clinsig_list)):
                        s = curr_clinsig_list[i]
                        if s != "":
                            try:
                                s_af = bp_af_adj_dict[d["bp_list"][i]]
                            except:
                                # The major allele was replaced, no score available for the correct substitution
                                continue
                            clinsig_list.append(s)
                            clinsig_af.append(s_af)

                    # Saving indices of syn and non-syn bps
                    syn_idx = []
                    nonsyn_idx = []
                    for i in range(len(d["bp_list"])):
                        ref_aa = d["aa_ref"]
                        alt_bp = d["bp_list"][i]
                        alt_aa = codon_table[alt_bp.upper()]
                        if alt_aa == ref_aa:
                            syn_idx.append(i)
                        else:
                            nonsyn_idx.append(i)

                    # Summing the AC per population
                    for i in range(len(ac_str)):
                        ac = ac_str[i]
                        ac_sum[i] += sum(d[ac])
                        # Summing syn and non-syn separately
                        ac_sum_syn[i] += sum(np.array(d[ac])[syn_idx])
                        ac_sum_nonsyn[i] += sum(np.array(d[ac])[nonsyn_idx])

                    # Averaging the AN per population, to do that, gathering all an to a list
                    for i in range(len(an_str)):
                        an = an_str[i]
                        (an_list[i]).extend(d[an])

                    # Averaging the MAF per population, to do that: gathering all maf!=0 to a list
                    for i in range(len(an_str)):
                        ac = ac_str[i]
                        an = an_str[i]
                        for j in range(len(d[ac])):
                            if d[an][j] != 0:
                                pop_maf = d[ac][j] / float(d[an][j])
                                if pop_maf != 0:
                                    if j in syn_idx:
                                        pop_maf_syn_list[i].append(pop_maf)
                                    else:
                                        pop_maf_nonsyn_list[i].append(pop_maf)
                                    pop_maf_list[i].append(pop_maf)

            features_dict['state_id'] = state_id
            features_dict['domain_name'] = domain

            for i in range(len(AMINO_ACIDS)):
                if aa_alt_prob[i] > 0:
                    aa_alt_prob_avg[i] = aa_alt_prob[i] / float(aa_alt_hist[i])

            features = [

                ExAC_MAF_features(sites_aa_num, sites_aa_alter_num, maf_list),
                ExAC_population_features(pop_maf_list, pop_maf_syn_list, pop_maf_nonsyn_list),
                ExAC_count_features(sites_aa_num, sites_aa_alter_num, sites_snp_num, sites_snp_alter_num),
                ExAC_rareSNP_features(sites_snp_alter_num, rare_5_num, rare_05_num, rare_005_num),

                conservation_features(phastCons_dict, phyloP_dict),

                sub_matrix_features(blosum62_list, weigted_blosum62_list, 'blosum'),
                sub_matrix_features(pam40_list, weigted_pam40_list, 'pam'),

                pseudo_dNdS_features(ref_seq, Nd, Sd),
                pfam_emission_prob_features(hmm_prob_dict, state),
                pfam_conserved_state_feature(state, con_states_dict),

                SIFT_features(sift_scores_list, weighted_sift_scores_list),
                PolyPhen_features(polyphen_scores_list, polyphen_pred_list, weighted_polyphen_scores_list),
                ClinVar_scores(clinsig_list, clinsig_af),

                entropy_features(maf_list),
                instance_individuals_100way_change_features(maf_list, aa_ref_hist, jsd100way_list),

                aa_identity_features(aa_ref_hist, "aa_ref"),
                major_allele_charge(aa_ref_hist),
                major_allele_hydrophobicity(aa_ref_hist),
                major_allele_volume(aa_ref_hist),
                major_allele_functional_group(aa_ref_hist),
                major_allele_propensity(aa_ref_hist),
                major_allele_h_bonds(aa_ref_hist),

                aa_identity_features(aa_alt_hist, "aa_alt_cnt"),
                aa_identity_features(aa_alt_prob_avg, "aa_alt_avg_freq"),

                sub_diff_hydrophobicity(ref_alt_pairs),
                sub_diff_volume(ref_alt_pairs),
                sub_diff_functional_group(ref_alt_pairs),
                sub_diff_propensity(ref_alt_pairs),
                sub_diff_h_bonds(ref_alt_pairs),

                spider_solvent_acc_pred(spider_dict),
                spider_contact_number_pred(spider_dict),
                spider_angles_pred(spider_dict),
                spider_struct_pred(spider_dict),
                spider_half_sphere_exposure_pred(spider_dict),

                whole_domain_conservation(states_dict),
                domain_location_features(state, max([k for k in states_dict if isinstance(k, int)])),
                protein_location_features(protein_pos_list, protein_len_list)

            ]

            for feature in features:
                features_dict.update(feature)
            features_list.append(features_dict)

    domains_features_df = pd.DataFrame(features_list)
    domains_features_df = domains_features_df.set_index('state_id')
    domains_features_df = domains_features_df.sort_index()

    to_lowercase_fields = 'maf_nonsyn_AMR maf_syn_SAS maf_NFE maf_nonsyn_EAS maf_AMR maf_nonsyn_FIN maf_syn_AFR maf_syn_EAS maf_syn_OTH maf_nonsyn_SAS maf_nonsyn_AFR maf_syn_FIN maf_FIN maf_nonsyn_NFE maf_EAS maf_AFR maf_nonsyn_OTH maf_SAS maf_OTH maf_syn_AMR maf_syn_NFE'.split()
    domains_features_df = domains_features_df.rename(columns={x: x.lower() for x in to_lowercase_fields})
    domains_features_df.to_csv(OUTPUT_CSV, sep='\t', index_label='')

