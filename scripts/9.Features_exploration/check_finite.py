import pandas as pd
pd.options.mode.use_inf_as_na = True

TSV = 'positions_features_fixed.csv'
TSV_ANAT = '/home/vineetb/git_checkouts/MyResearch/ExAC/9.Features_exploration/ig.csv'


if __name__ == '__main__':

    df = pd.read_csv(TSV, sep='\t', index_col=0)
    print(df.shape)

    for c in df.columns:
        na_series = df[c].isna()

        if any(na_series):
            na_indices = na_series[na_series].index.values
            print(f'column {c} - {len(na_indices)} na - {na_indices}')

    df_anat = pd.read_csv(TSV_ANAT, sep='\t', index_col=0)
    df_anat = df_anat[df_anat.domain_name=='ig']
    print(df_anat.shape)

    extra_cols = [c for c in df.columns if c not in df_anat.columns]
    print('extra cols = ' + str(extra_cols))

    missing_cols = [c for c in df_anat.columns if c not in df.columns]
    print('missing cols = ' + str(missing_cols))

    for c in df_anat.columns:
        na_series = df_anat[c].isna()

        if any(na_series):
            na_indices = na_series[na_series].index.values
            print(f'column {c} - {len(na_indices)} na - {na_indices}')
