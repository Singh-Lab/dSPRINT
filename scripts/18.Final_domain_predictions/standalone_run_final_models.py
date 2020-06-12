import os.path
import pandas as pd
import pickle
import dsprint.models

try:
    snakemake
except NameError:
    import sys
    if len(sys.argv) != 3:
        print('Usage: <script> <input_csv> <output_csv>')
        sys.exit(0)

    INPUT_CSV, OUTPUT_CSV = sys.argv[1:]
else:
    INPUT_CSV = snakemake.input.input_csv
    OUTPUT_CSV = snakemake.output.output_csv


MODELS_REQ_SCALING = ['SVM', 'KNN', 'Logistic', 'NN']
models_dir = os.path.dirname(dsprint.models.__file__)

with open(os.path.join(models_dir, 'level1', 'scaler.pik'), 'rb') as f:
    SCALER = pickle.load(f, encoding='latin1')


def predict(name, pik_model, classifier_method, data):
    if classifier_method in MODELS_REQ_SCALING:
        data = data.copy()
        data[:] = SCALER.transform(data)

    prob = pik_model.predict_proba(data)
    if classifier_method != 'NN':
        prob = prob[:, 1].tolist()

    return pd.Series(prob, name=name, index=data.index)


if __name__ == '__main__':

    MODELS = ['XGB', 'RF', 'SVM', 'Logistic', 'NN']
    LIGANDS = ['dna', 'rna', 'ion', 'peptide', 'sm']

    features = pd.read_csv(INPUT_CSV, sep='\t', index_col=0)
    features.drop('domain_name', axis=1, inplace=True)

    # ------------
    # LAYER 1
    # ------------
    all_predictions = []
    for model in MODELS:
        for ligand in LIGANDS:
            with open(os.path.join(models_dir, 'level1', f'{ligand}_{model}.pik'), 'rb') as f:
                pik_model = pickle.load(f, encoding='latin1')
            predictions = predict(f'{model}_{ligand}_prob', pik_model, model, features)
            all_predictions.append(predictions)

    predictions_level1 = pd.concat(all_predictions, axis=1)

    # ------------
    # LAYER 2
    # ------------
    # Mapping from 2nd_layer ligand type to a 2-tuple: (<1st_layer_ligand_types>, <1st_layer_models>)
    ligand_to_1st_level = {
        'sm':      (LIGANDS, MODELS),
        'dna':     (['dna'], MODELS),
        'ion':     (LIGANDS, ['XGB']),
        'peptide': (LIGANDS, ['XGB']),
        'rna':     (LIGANDS, ['XGB'])
    }

    all_predictions = []
    for ligand in LIGANDS:
        with open(os.path.join(models_dir, 'level2', f'{ligand}_XGB.pik'), 'rb') as f:
            pik_model = pickle.load(f, encoding='latin1')

        # All feature names we wish to extract from the 1st level predictions
        ligands_level1, models_level1 = ligand_to_1st_level[ligand]
        feature_names_level1 = [f'{model}_{ligand}_prob' for ligand in ligands_level1 for model in models_level1]

        # Get level 1 predictions, create a combined dataframe with features + level 1 predictions
        features_level1 = predictions_level1[feature_names_level1]
        features_level1 = features_level1.merge(features, left_index=True, right_index=True)

        predictions = predict(ligand, pik_model, 'XGB', features_level1)
        all_predictions.append(predictions)

    # combine all predictions into a dataframe, where every column is the ligand name
    df = pd.concat(all_predictions, axis=1)

    """
    Since our model predictions are only based on features present as columns, we have thus far ignored the domain name
    and the the 1-indexed state, by moving them to the index. Now we move them back in as 'domain', 'match_state'.
    We also add in the 'ligand_type' and 'binding_score' columns by unpivoting columns to rows.
    The resulting dataframe has domain/match_state/ligand_type/binding_score, suitable for ingestion at 
    protdomain.princeton.edu
    """

    df = pd.melt(df.reset_index(), id_vars='index', var_name='ligand_type', value_name='binding_score')
    df['domain'], df['match_state'] = df['index'].str.split('_').str
    df = df.drop('index', axis=1).astype({'match_state': int})
    df.to_csv(OUTPUT_CSV, sep='\t')
