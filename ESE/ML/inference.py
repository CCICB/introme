from itertools import chain
import os, pickle, sys

from vcf2pandas import vcf2pandas
import pandas as pd
import numpy as np

INFO_FIELDS = {
    # ./-/+ --> 0:. 1:-/+
    'AG_Created': 'AGcheck_AG_Created',
    'AG_Lost': 'AGcheck_AG_Lost',
    'GT_Created': 'AGcheck_GT_Created',
    'GT_Lost': 'AGcheck_GT_Lost',
    # SNV/INDEL/INSDEL --> one hot encode
    'Variant_Type': 'AGcheck_Variant_Type',

    'Branchpointer_Prob': 'Branchpointer_Prob',
    'Branchpointer_U2_Binding_Energy': 'Branchpointer_U2_Binding_Energy',
    'Branchpointer_max_Prob': 'Branchpointer_max_Prob',
    'Branchpointer_max_U2_Binding_Energy': 'Branchpointer_max_U2_Binding_Energy',
    'Branchpointer_options': 'Branchpointer_options',

    'ESE_A1_Hazeem_alt': 'ESE_A1_Hazeem_alt',
    'ESE_A1_Hazeem_diff': 'ESE_A1_Hazeem_diff',
    'ESE_A1_Hazeem_july_alt': 'ESE_A1_Hazeem_july_alt',
    'ESE_A1_Hazeem_july_diff': 'ESE_A1_Hazeem_july_diff',
    'ESE_A1_neuBG_alt': 'ESE_A1_neuBG_alt',
    'ESE_A1_neuBG_diff': 'ESE_A1_neuBG_diff',
    'ESE_A1_winBG_alt': 'ESE_A1_winBG_alt',
    'ESE_A1_winBG_diff': 'ESE_A1_winBG_diff',
    'ESE_SRSF1_alt': 'ESE_SRSF1_alt',
    'ESE_SRSF1_diff': 'ESE_SRSF1_diff',
    'ESE_SRSF1_igM_alt': 'ESE_SRSF1_igM_alt',
    'ESE_SRSF1_igM_diff': 'ESE_SRSF1_igM_diff',
    'ESE_SRSF2_alt': 'ESE_SRSF2_alt',
    'ESE_SRSF2_diff': 'ESE_SRSF2_diff',
    'ESE_SRSF5_alt': 'ESE_SRSF5_alt',
    'ESE_SRSF5_diff': 'ESE_SRSF5_diff',
    'ESE_SRSF6_alt': 'ESE_SRSF6_alt',
    'ESE_SRSF6_diff': 'ESE_SRSF6_diff',

    'MMSplice_alt_acceptor': 'MMSplice_alt_acceptor',
    'MMSplice_alt_acceptor_intron': 'MMSplice_alt_acceptor_intron',
    'MMSplice_alt_donor': 'MMSplice_alt_donor',
    'MMSplice_alt_donor_intron': 'MMSplice_alt_donor_intron',
    'MMSplice_alt_exon': 'MMSplice_alt_exon',
    'MMSplice_delta_logit_PSI': 'MMSplice_delta_logit_PSI',
    'MMSplice_pathogenicity': 'MMSplice_pathogenicity',
    'MMSplice_ref_acceptor': 'MMSplice_ref_acceptor',
    'MMSplice_ref_acceptor_intron': 'MMSplice_ref_acceptor_intron',
    'MMSplice_ref_donor': 'MMSplice_ref_donor',
    'MMSplice_ref_donor_intron': 'MMSplice_ref_donor_intron',
    'MMSplice_ref_exon': 'MMSplice_ref_exon',

    'Pangolin_Gain': 'Pangolin_Gain',
    'Pangolin_Loss': 'Pangolin_Loss',

    'SPIP_Range_Max': 'SPIP_Range_Max',
    'SPIP_Range_Min': 'SPIP_Range_Min',
    'SPIP_Score': 'SPIP_Score',
    'SPiCE_Prob': 'SPIP_SPiCE_Prob',

    'DP_AG': 'SpliceAI_DP_AG',
    'DP_AL': 'SpliceAI_DP_AL',
    'DP_DG': 'SpliceAI_DP_DG',
    'DP_DL': 'SpliceAI_DP_DL',
    'DS_AG': 'SpliceAI_DS_AG',
    'DS_AG_ALT': 'SpliceAI_DS_AG_ALT',
    'DS_AL': 'SpliceAI_DS_AL',
    'DS_AL_ALT': 'SpliceAI_DS_AL_ALT',
    'DS_DG': 'SpliceAI_DS_DG',
    'DS_DG_ALT': 'SpliceAI_DS_DG_ALT',
    'DS_DL': 'SpliceAI_DS_DL',
    'DS_DL_ALT': 'SpliceAI_DS_DL_ALT',

    'AccGainP': 'Spliceogen_AccGainP',
    'AccLossP': 'Spliceogen_AccLossP',
    'DonGainP': 'Spliceogen_DonGainP',
    'DonLossP': 'Spliceogen_DonLossP',
    'ESEmaxAlt': 'Spliceogen_ESEmaxAlt',
    'ESEmaxRef': 'Spliceogen_ESEmaxRef',
    'ESSminAlt': 'Spliceogen_ESSminAlt',
    'ESSminRef': 'Spliceogen_ESSminRef',
    'mesAccAlt': 'Spliceogen_mesAccAlt',
    'mesAccRef': 'Spliceogen_mesAccRef',
    'mesDonAlt': 'Spliceogen_mesDonAlt',
    'mesDonRef': 'Spliceogen_mesDonRef',

    # U12/.      --> 0:. 1:U12
    'Intron_Type': 'U12_Intron_Type',
    # score/dots --> 0:. else score
    'U12_score': 'U12_score',
    # ./-/+      --> 0:. 1:-/+
    # 'U12_strand': 'U12_strand',

    # [('acceptor_canonical',) ('donor_region',) '.' ('donor_exonic',)
    #  ('acceptor_exonic',) ('donor_canonical',) ('branchpoint_region',)
    #  ('acceptor_region',)]
    'gene_regions': 'GENEINFO_gene_regions',

    ### Gencode file
    # intronic/exonic
    'Gene_Location': 'GENEINFO_is_intronic',
    # Exclude
    # 'gene': 'gene',

    # Exclude
    # [('protein_coding',) ('protein_coding', 'lncRNA')
    #  ('lncRNA', 'protein_coding')
    #  ('lncRNA', 'transcribed_unprocessed_pseudogene') ('lncRNA',) etc.
    # 'gene_type': 'gene_type',

    # [('-',) ('+',) ('+', '-') ('-', '+')]
    'strand': 'GENEINFO_strand',
    }

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


def main(model, df, columns, outfile):
    print(df.head(2))
    df_original = df.copy()
    
    ese        = [col for col in df.columns if col.startswith('INFO:ESE_')]
    mmsplice   = [col for col in df.columns if col.startswith('INFO:MMSplice_')]
    pangolin   = [col for col in df.columns if col.startswith('INFO:Pangolin_')]
    spip       = [col for col in df.columns if col.startswith('INFO:SPIP_')]
    spliceai   = [col for col in df.columns if col.startswith('INFO:SpliceAI_DS')]
    spliceogen = [col for col in df.columns if col.startswith('INFO:Spliceogen_')]

    ensemble_scores = [ese, mmsplice, pangolin, spliceai, spliceogen, spip]

    ## 1. Preprocessing data in the same way as training... next time use sklearn pipelines
    DUMMY_COLS = ['INFO:AGcheck_Variant_Type', 'INFO:GENEINFO_gene_regions']
    GENE_REGIONS = [f"INFO:GENEINFO_gene_regions_{suffix}" for suffix in [
                        ".",
                        "('acceptor_canonical',)",
                        "('acceptor_exonic',)",
                        "('acceptor_region',)",
                        "('branchpoint_region',)",
                        "('donor_canonical',)",
                        "('donor_exonic',)",
                        "('donor_region',)"
                    ]]
    VARIANT_TYPE = [f"INFO:AGcheck_Variant_Type_{suffix}" for suffix in [
                        "('INDEL',)",
                        "('INSDEL',)",
                        "('SNV',)",
                    ]]


    for col in [col for col in df.columns if col.startswith('INFO:AGcheck_')]:
        if not col.endswith('Created') and not col.endswith('Lost'):
            continue
        # Convert AGcheck_lost/created to 0 if dot, and 1 if minus OR plus strand
        df[col] = (df[col] != '.')

    df['INFO:U12_Intron_Type'] = (df['INFO:U12_Intron_Type'] == 'U12')
    df['INFO:GENEINFO_is_intronic'] = (df['INFO:GENEINFO_is_intronic'] == ('intronic',))

    df['INFO:GENEINFO_strand_minus'] = df['INFO:GENEINFO_strand'].apply(lambda x: '-' in x)
    df['INFO:GENEINFO_strand_plus'] = df['INFO:GENEINFO_strand'].apply(lambda x: '+' in x)
    df = df.drop(columns=['INFO:GENEINFO_strand'])
    
    # 1a. Get dummies for AGcheck variant type AND GENEINFO gene regions
    dummys = pd.get_dummies(df[DUMMY_COLS], columns=DUMMY_COLS).reindex(
        columns=GENE_REGIONS+VARIANT_TYPE, fill_value=False
    )

    df = pd.concat([dummys, df.drop(columns=DUMMY_COLS)], axis=1)
    print(len(df.columns))

    print(df.head(2))

    print(set(df.columns) - set(columns))
    print(set(columns) - set(df.columns))

    df = df.reindex(columns=columns)
    print(df.columns)
    print(len(df.columns))

    # agcheck    = [col for col in df.columns if col.startswith('INFO:AGcheck_')]
    # bpter      = [col for col in df.columns if col.startswith('INFO:Branchpointer_')]
    # u12        = [col for col in df.columns if col.startswith('INFO:U12_')]
    # geneinfo   = [col for col in df.columns if col.startswith('INFO:GENEINFO_')]
    # info       = [agcheck, bpter, u12, geneinfo]


    # # 1b. Sanitise further
    # scores = list(chain.from_iterable(ensemble_scores))
    # scores.extend(list(chain.from_iterable(info)))

    print(df.head(2))
    features = assert_and_convert_single_float_tuples_allow_dot(df,
                                                                dots_largeminus=['INFO:SpliceAI_DP',
                                                                                'INFO:SPIP'],
                                                                dots_zeros=['INFO:U12',
                                                                            'INFO:Branchpointer',
                                                                            'INFO:SpliceAI_DS',
                                                                            'INFO:Spliceogen'],
                                                                nans_zeros=[],
                                                                nans_large_minus=['INFO:SPIP'])
    print("before:", len(features))

    bad_features = features[features.isin(["."]).any(axis=1)]
    features = features[~features.isin(["."]).any(axis=1)]
    print("after:", len(features))

    # print("Removed bad features:")
    # print(bad_features.to_string())

    y_scores = model.predict_proba(features)[:, 1]
    scores_ser = pd.Series(y_scores, index=features.index).astype(object)
    scores_full = scores_ser.reindex(df_original.index, fill_value=None)

    print(scores_full.head())

    df_original['introme_score'] = scores_full
    df_original.to_csv(outfile, sep="\t")


if __name__ == "__main__":
    import json

    USAGE = f"""
    Usage: {sys.argv[0]} model_path.pkl columns.json splicing_anno.vcf output_path.tsv
    """
    if len(sys.argv) != 5:
        raise ValueError(USAGE)

    model_path = sys.argv[1]
    columns_path = sys.argv[2]
    splicing_anno_vcf_path = sys.argv[3]
    output_path = sys.argv[4]

    with open(model_path, 'rb') as file:
        model = pickle.load(file)

    with open(columns_path, 'rb') as file:
        columns = json.load(file)

    print(columns, len(columns), type(columns))
    # exit(0)

    df = vcf2pandas(splicing_anno_vcf_path,
                    remove_empty_columns=False,
                    info_fields=INFO_FIELDS)

    with open(output_path, 'w') as outfile:
        main(model, df, columns, outfile)
