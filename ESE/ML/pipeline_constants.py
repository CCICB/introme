from typing import Any, TypeAlias

from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier, GradientBoostingClassifier

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

###############################################################################
# Training Settings below

ENSEMBLE_SCORE_COLS = {
    "ese": "INFO:ESE_",
    "mmsplice": "INFO:MMSplice_",
    "pangolin": "INFO:Pangolin_",
    "spip": "INFO:SPIP_",
    "spliceai": "INFO:SpliceAI_DS",
    "spliceogen": "INFO:Spliceogen_",
}

# These are the main scores from each tool, usually splice change prob or log prob.
# Apply max of abs to all specified columns
RAW_SCORE_COLS = {
    "spliceai": ['INFO:SpliceAI_DS_AG', 'INFO:SpliceAI_DS_AL', 'INFO:SpliceAI_DS_DG', 'INFO:SpliceAI_DS_DL'],
    "spip": ['INFO:SPIP_SPiCE_Prob'],
    "pangolin": ['INFO:Pangolin_Gain', 'INFO:Pangolin_Loss'],
    "spliceogen": ['INFO:Spliceogen_AccGainP', 'INFO:Spliceogen_AccLossP', 'INFO:Spliceogen_DonGainP', 'INFO:Spliceogen_DonLossP'],
    # "The main score is predicted by MMSplice, which shows the effect of the variant on the inclusion level (PSI percent spliced in) of the exon.
    # If delta_logit_psi is bigger than 2 or smaller than -2, the effect of variant can be considered strong."
    "mmsplice": ['INFO:MMSplice_delta_logit_PSI'],
}

ClassifierCtor: TypeAlias = type[RandomForestClassifier | HistGradientBoostingClassifier | GradientBoostingClassifier]
ClassifierSpec: TypeAlias = tuple[ClassifierCtor, dict[str, Any]]

CLASSIFIERS: dict[str, ClassifierSpec] = {
    "RandomForest": (RandomForestClassifier, {'n_jobs': -1, 'min_samples_leaf': 5}),
    "XGBoost": (GradientBoostingClassifier, {'min_samples_leaf': 5}),
    "HistGradientBoosting": (HistGradientBoostingClassifier, {'min_samples_leaf': 5}),
}