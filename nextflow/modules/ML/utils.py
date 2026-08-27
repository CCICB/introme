import pandas as pd
import numpy as np

from pipeline_constants import (
    DUMMY_COLS,
    GENE_REGIONS,
    VARIANT_TYPE,
)

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
                                                     dots_largeminus: list[str] = [],
                                                     dots_zeros: list[str] = [],
                                                     nans_zeros: list[str] = [],
                                                     nans_large_minus: list[str] = []) -> pd.DataFrame:
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
                # return (len(x) == 1 and (isinstance(x[0], float) or x[0] is None))
                return ((len(x) == 1) and (isinstance(x[0], float) or x[0] is None)) or (len(x) > 1 and all(elem == x[0] for elem in x[1:]))

            # If any row fails the validity check, raise an error
            valid_series = df[col].apply(is_valid_tuple)
            if not valid_series.all():
                # Extract offending values
                offending_values = df[col][~valid_series].unique()
                raise ValueError(
                    f"Column '{col}' contains invalid tuple values: {offending_values}"
                )

        # Convert single-float tuples to that float
        def convert_tuple_to_singleton(x):
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

        df[col] = df[col].apply(convert_tuple_to_singleton)
        df[col] = pd.to_numeric(df[col], errors='coerce')

    return df

def filter_bad_features(X_split, y_split, feature_cols):
    """
    Remove rows where any feature is "." (indicating unscorable by that tool)
    """
    X_subset = X_split[feature_cols].copy()
    initial_len = len(X_subset)

    bad_feature_mask = X_subset.isna().any(axis=1) | (X_subset == '.').any(axis=1)
    X_clean, y_clean = X_subset[~bad_feature_mask], y_split[~bad_feature_mask]

    dropped_len = initial_len - len(X_clean)
    if dropped_len > 0:
        print(f"cleaned rows. before: {initial_len} | after: {len(X_clean)} | dropped: {dropped_len}")
    return X_clean, y_clean

def fmt_float32(value, dp: int=-1):
    """Format a number like vcfgo's fmtFloat32."""
    if value == "." or value is None:
        return value

    # Preserve integer INFO fields exactly rather than converting them to float32.
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return str(value)

    v = np.float32(value)

    if dp >= 0:
        val = f"{v:.{dp}f}"
    else:
        if v > 0.02 or v < -0.02:
            val = f"{v:.4f}"
        else:
            val = f"{v:.5g}"

    val = val.rstrip("0").rstrip(".")

    if val in ("", "-"):
        val = "0"

    return val

def is_numeric_info_column(series):
    def is_numeric_value(x):
        if x is None or x == ".":
            return True

        if isinstance(x, (float, np.floating)) and np.isnan(x):
            return True

        return isinstance(x, (int, float, np.integer, np.floating)) \
            and not isinstance(x, bool)

    return series.map(is_numeric_value).all()
