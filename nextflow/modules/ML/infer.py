from pathlib import Path
import json
import pickle

import pandas as pd
import numpy as np

from utils import preprocess, filter_bad_features

def infer_main(model_path: Path, df: pd.DataFrame, columns_path: Path, output_path: Path):
    """
    Preprocess according to inference mode, run inference on supplied model and columns
    and save output TSV.

    1-5. Preprocess
    6. (GENERATED) Reorder columns to match training order
        - Given by param (combinations) for inference, should be the training columns of used model.
    """
    with open(model_path, 'rb') as file:
        model = pickle.load(file)

    with open(columns_path, 'r') as file:
        columns = json.load(file)

    df_original = df.copy()
    df = preprocess(train_mode=False, df=df)

    print(f"vcf has {len(df.columns)} columns (including dummies)")

    print(df.head(2))

    print("excluded columns:", set(df.columns) - set(columns))
    print("expected but unpresent columns:", set(columns) - set(df.columns))

    # features = df[columns]

    features, _ = filter_bad_features(df, pd.Series([0]*len(df)), columns)

    if len(features) == 0:
        print("No valid rows after filtering for inference. Saving empty output.")
        df['introme_score'] = np.nan
        df.to_csv(output_path, sep='\t', index=False)
        return

    y_scores = model.predict_proba(features)[:, 1]
    scores_series = pd.Series(y_scores, index=features.index).astype(object)
    scores_full = scores_series.reindex(df_original.index, fill_value=None)

    # for col starting with "INFO:" if tuple but only length one, unpack to just that value
    for col in df_original.columns:
        if col.startswith("INFO:") and df_original[col].apply(lambda x: (isinstance(x, tuple) and len(x) == 1) or x is None or x == ".").all():
            df_original[col] = df_original[col].apply(lambda x: x[0] if isinstance(x, tuple) and len(x) == 1 else x)
    df_original['introme_score'] = scores_full
    df_original.to_csv(output_path, sep='\t', index=False)
