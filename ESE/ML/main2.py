from itertools import chain
import json
import os, pickle, sys
from typing import Literal

from vcf2pandas import vcf2pandas
import pandas as pd
import numpy as np

from pipeline_constants import (
    DUMMY_COLS,
    GENE_REGIONS,
    INFO_FIELDS,
    VARIANT_TYPE,
)

"""
Functions:
- preprocess()
    - assert_and_convert_single_float_tuples_allow_dot()
- subset_features()
- train_main()
- infer_main()
"""

def preprocess(train_mode: bool, df: pd.DataFrame) -> pd.DataFrame:
    """
    Training mode:
    1. Remove conflicting class and dummify the other class columns
    2. Dummyify AGcheck_Variant_Type and GENEINFO_gene_regions, then drop original columns
    3. (Common) AG_check lost/created dots -> 0/1's
    4. (Common) U12 intronic & Strand dots -> 0/1's
    5. (Common) dots/NaNs in numeric columns -> 0's or large negative values (configurable)
    
    Inference mode:
    1. (SKIP because there are no ground truth labels)
    2. (CONSTANT) Generate the seen dummy columns from pipeline_constants
    3. (Common) AG_check lost/created dots -> 0/1's
    4. (Common) U12 intronic & Strand dots -> 0/1's
    5. (Common) dots/NaNs in numeric columns -> 0's or large negative values (configurable)
    """
    if train_mode:
        # 1. Remove conflicting class and dummify the other class columns
        # (e.g. if label is binary, drop one of the two columns and rename the other to "label")
        # 2. Dummyify AGcheck_Variant_Type and GENEINFO_gene_regions, then drop original columns
        df = df[df['ID'] != 'SVDVConflicting']
        df = pd.get_dummies(df, columns=['ID'])

        df = pd.get_dummies(df, columns=['INFO:AGcheck_Variant_Type'])
        df = pd.get_dummies(df, columns=['INFO:GENEINFO_gene_regions'])
    else:
        # 1. (SKIP because there are no ground truth labels)
        # 2. (CONSTANT) Generate the seen dummy columns from pipeline_constants
        dummys = pd.get_dummies(df[DUMMY_COLS], columns=DUMMY_COLS).reindex(
            columns=GENE_REGIONS + VARIANT_TYPE,
            fill_value=False,
        )
        df = pd.concat([dummys, df.drop(columns=DUMMY_COLS)], axis=1)

    # 3. (Common) AG_check lost/created dots -> 0/1's
    for col in [col for col in df.columns if col.startswith('INFO:AGcheck_')]:
        if not col.endswith('Created') and not col.endswith('Lost'):
            continue
        # Convert AGcheck_lost/created to 0 if dot, and 1 if minus OR plus strand
        df[col] = (df[col] != '.')

    # 4. (Common) U12 intronic & Strand dots -> 0/1's
    df['INFO:U12_Intron_Type'] = (df['INFO:U12_Intron_Type'] == 'U12')
    df['INFO:GENEINFO_is_intronic'] = (df['INFO:GENEINFO_is_intronic'] == ('intronic',))

    df['INFO:GENEINFO_strand_minus'] = df['INFO:GENEINFO_strand'].apply(lambda x: '-' in x)
    df['INFO:GENEINFO_strand_plus'] = df['INFO:GENEINFO_strand'].apply(lambda x: '+' in x)
    df = df.drop(columns=['INFO:GENEINFO_strand'])

    # 5. (Common) dots/NaNs in numeric columns -> 0's or large negative values (configurable)
    features = assert_and_convert_single_float_tuples_allow_dot(df,
                                                                dots_largeminus=['INFO:SpliceAI_DP',
                                                                                'INFO:SPIP'],
                                                                dots_zeros=['INFO:U12',
                                                                            'INFO:Branchpointer',
                                                                            'INFO:SpliceAI_DS',
                                                                            'INFO:Spliceogen',
                                                                            'INFO:MMSplice',
                                                                            'INFO:Pangolin'],
                                                                nans_zeros=[],
                                                                nans_large_minus=['INFO:SPIP'])

    return features

def assert_and_convert_single_float_tuples_allow_dot(df: pd.DataFrame, *,
                                                     dots_largeminus: list[str],
                                                     dots_zeros: list[str],
                                                     nans_zeros: list[str],
                                                     nans_large_minus: list[str]) -> pd.DataFrame:
    """
    For each column in df:
      - If the column contains any tuples, assert that every tuple is either:
          * length-1 and contains a float, or
          * is None, or
          * equals the string "."
      - Convert such valid single-float tuples to the float value.
      - Leave the "." string as-is instead of trying to convert it to float.
    Returns a copy of df with the converted columns (dtype may remain 'object' 
    if "." strings are present).
    """
    df = df.copy()
    
    for col in df.columns:
        if not col.startswith("INFO:"):
            continue
        else:
            print(f"Processing column: {col}")
        # Check if this column contains *any* tuples
        has_tuple = df[col].apply(lambda x: isinstance(x, tuple)).any()
        
        if has_tuple:
            # Define a helper to check if each value is valid
            def is_valid_tuple(x):
                # Valid if it's None, or the special ".", or a single-float tuple
                if x is None or x == "." or not isinstance(x, tuple):
                    return True
                if (x[0] is None):
                    print(col, x[0])
                return (len(x) == 1 and (isinstance(x[0], float) or x[0] is None))
            
            # If any row fails the validity check, raise an error
            valid_series = df[col].apply(is_valid_tuple)
            if not valid_series.all():
                # Extract offending values
                offending_values = df[col][~valid_series].unique()
                raise ValueError(
                    f"Column '{col}' contains invalid tuple values: {offending_values}"
                )
            
        # Convert single-float tuples to that float
        def convert_tuple_to_float(x):
            if isinstance(x, tuple):
                return sanitise(x[0])  # we've already asserted it's length-1 float
            return sanitise(x)
        
        def sanitise(x):
            if not (str(x) in ["", ".", "NaN", "nan", "None"]):
                return x
            # because spliceai can't do multinucleotide to multinucleotide variants.
            if any(col.startswith(x) for x in dots_largeminus) and (str(x) == "." or str(x) == ""):
                return -99999
            if any(col.startswith(x) for x in dots_zeros) and (str(x) == "." or str(x) == ""):
                return 0
            if any(col.startswith(x) for x in nans_zeros) and str(x) in ["NaN", "nan", "None"]:
                return 0
            if any(col.startswith(x) for x in nans_large_minus) and str(x) in ["NaN", "nan", "None"]:
                return -99999
            return x
        
        df[col] = df[col].apply(convert_tuple_to_float)
    
    return df

# def subset_features(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
#     """
#     Group df columns by prediction tools, also into 'ensemble' and 'info' (auxiliary info)
#     This is a given param (list, and category name mapping)
    
#     Train:

#     7. (COMMON) Return df with selected columns, and the list of selected columns
        
#     Inference:

#     7. (COMMON) Return df with selected columns, and the list of selected columns
#     """
#     df = df.reindex(columns=columns)
#     return df

def make_col_combinations(df: pd.DataFrame) -> dict[str, list[str]]:
    """
    Generate feature combinations based on the DataFrame columns.
    """
    ensemble_scores = [
        ("ese",         [c for c in df.columns if c.startswith("INFO:ESE_")]),
        ("mmsplice",    [c for c in df.columns if c.startswith("INFO:MMSplice_")]),
        ("pangolin",    [c for c in df.columns if c.startswith("INFO:Pangolin_")]),
        ("spip",        [c for c in df.columns if c.startswith("INFO:SPIP_")]),
        ("spliceai",    [c for c in df.columns if c.startswith("INFO:SpliceAI_DS")]),
        ("spliceogen",  [c for c in df.columns if c.startswith("INFO:Spliceogen_")]),
    ]

    ensemble_cols = [cols for _, cols in ensemble_scores]

    agcheck    = [col for col in df.columns if col.startswith('INFO:AGcheck_')]
    bpter      = [col for col in df.columns if col.startswith('INFO:Branchpointer_')]
    u12        = [col for col in df.columns if col.startswith('INFO:U12_')]
    geneinfo   = [col for col in df.columns if col.startswith('INFO:GENEINFO_')]
    info       = [agcheck, bpter, u12, geneinfo]

    combinations: dict[str, list[str]] = {
        "all_info": list(chain.from_iterable(ensemble_cols)) + list(chain.from_iterable(info)),
        "all_noinfo": list(chain.from_iterable(ensemble_cols)),
        "notools_info": list(chain.from_iterable(info)),
    }

    for name, cols in ensemble_scores:
        combinations[name + "_info"] = cols + list(chain.from_iterable(info))
        combinations[name + "_noinfo"] = cols

    return combinations

def train_main(model, df, columns, outfile):
    """
    Preprocess according to train mode, train various models and save the best model weight/model/columns.
    1-5. Preprocess
    6. Subset to columns used for training

    For each model, for example:
        Model 1: Random Forest Classifier
        Model 2: Logistic Regression
        Model 3: XGBoost Classifier
        Model 4: HGBoost Classifier
    Do the following:
        a) Hyperparameter tuning
        b) Different random seeds
        c) Different feature subset (e.g. SpliceAI only, Spliceogen only, etc.)
        d) Evaluate on validation set, save F1, AUPRC, etc.
    """
    df = preprocess(train_mode=True, df=df)
    ensemble_scores = [
        ("ese",         [c for c in df.columns if c.startswith("INFO:ESE_")]),
        ("mmsplice",    [c for c in df.columns if c.startswith("INFO:MMSplice_")]),
        ("pangolin",    [c for c in df.columns if c.startswith("INFO:Pangolin_")]),
        ("spip",        [c for c in df.columns if c.startswith("INFO:SPIP_")]),
        ("spliceai",    [c for c in df.columns if c.startswith("INFO:SpliceAI_DS")]),
        ("spliceogen",  [c for c in df.columns if c.startswith("INFO:Spliceogen_")]),
    ]

    ensemble_cols = [cols for _, cols in ensemble_scores]

    agcheck    = [col for col in df.columns if col.startswith('INFO:AGcheck_')]
    bpter      = [col for col in df.columns if col.startswith('INFO:Branchpointer_')]
    u12        = [col for col in df.columns if col.startswith('INFO:U12_')]
    geneinfo   = [col for col in df.columns if col.startswith('INFO:GENEINFO_')]
    info       = [agcheck, bpter, u12, geneinfo]

    # for l in ensemble_scores:
    #     print(len(l), l)

    # for l in info:
    #     print(len(l), l)
    
    # for col in df.columns:
    #     if col == "REF" or col == "ALT":
    #         continue
    #     print(f"Column: {col}")
    #     print("Unique values:", df[col].unique())
    #     print()

    # df[df['ID'] == 'SVDBConflicting']

    combinations: dict[str, list[str]] = {
        "all_info": list(chain.from_iterable(ensemble_cols)) + list(chain.from_iterable(info)),
        "all_noinfo": list(chain.from_iterable(ensemble_cols)),
        "notools_info": list(chain.from_iterable(info)),
    }

    for name, cols in ensemble_scores:
        combinations[name + "_info"] = cols + list(chain.from_iterable(info))

    for name, cols in combinations.items():
        print(f"Training with feature set: {name} ({len(cols)} features)")
        # cols = cols + ['ID_SVDBSplice']
        features = df[cols]
        target = df['ID_SVDBSplice']
        # features = features.drop(columns=['ID_SVDBSplice'])

        # print(f"{len(cols)=} {cols=}")

        # print("before:", len(features))
        bad_features = features[features.isin(["."]).any(axis=1)]
        features = features[~features.isin(["."]).any(axis=1)]
        # print("after:", len(features))

        # print("Removed bad features:")
        # print(bad_features.to_string())

        # ones = (target == 1).sum()
        # zeros = (target == 0).sum()

        # print(f'{target}, {ones=}, {zeros=}')

 

def infer_main(model, df, columns, outfile):
    """
    Preprocess according to inference mode, run inference on supplied model and columns
    and save output TSV.

    1-5. Preprocess
    6. (GENERATED) Reorder columns to match training order
        - Given by param (combinations) for inference, should be the training columns of used model.
    """
    df = preprocess(train_mode=False, df=df)

    print(f"vcf has {len(df.columns)} columns (including dummies)")

    print(df.head(2))

    print("excluded columns:", set(df.columns) - set(columns))
    print("expected but unpresent columns:", set(columns) - set(df.columns))

    features = df[columns]

    print("before:", len(features))

    bad_features = features[features.isin(["."]).any(axis=1)]
    features = features[~features.isin(["."]).any(axis=1)]
    print("after:", len(features))

    obj_cols = features.select_dtypes(include=["object"]).columns.tolist()
    print("Object dtype columns:", obj_cols)

    y_scores = model.predict_proba(features)[:, 1]
    scores_series = pd.Series(y_scores, index=features.index).astype(object)
    scores_full = scores_series.reindex(df.index, fill_value=None)

    df['introme_score'] = scores_full
    df.to_csv(outfile, sep='\t')

if __name__ == "__main__":
    USAGE = f"""
    Usage: {sys.argv[0]} [train|infer] model_path.pkl columns.json splicing_anno.vcf output_path.tsv
    """
    if len(sys.argv) != 6:
        raise ValueError(USAGE)

    is_train_mode = (sys.argv[1] == "train")
    model_path = sys.argv[2]
    columns_path = sys.argv[3]
    splicing_anno_vcf_path = sys.argv[4]
    output_path = sys.argv[5]

    with open(model_path, 'rb') as file:
        model = pickle.load(file)

    with open(columns_path, 'rb') as file:
        columns = json.load(file)

    print('columns loaded', columns, len(columns), type(columns))
    # exit(0)

    df = vcf2pandas(splicing_anno_vcf_path,
                    remove_empty_columns=is_train_mode,
                    info_fields=INFO_FIELDS)

    if is_train_mode:
        train_main(model, df, columns, output_path)
    else:
        with open(output_path, 'w') as outfile:
            infer_main(model, df, columns, outfile)