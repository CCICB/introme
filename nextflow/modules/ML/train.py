from itertools import chain
from pathlib import Path
import pickle
import json

import numpy as np
import pandas as pd

# Machine Learning
from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier, GradientBoostingClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score, precision_recall_curve, auc
from sklearn.model_selection import GroupKFold

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from pipeline_constants import (
    DUMMY_COLS,
    GENE_REGIONS,
    INFO_FIELDS,
    VARIANT_TYPE,
    ENSEMBLE_SCORE_COLS,
    RAW_SCORE_COLS,
    CLASSIFIERS,
)

from utils import preprocess, filter_bad_features

TARGET_COL = 'ID_SVDBSplice'
CHROM_COL = 'CHROM'
DEFAULT_TEST_CHROMS = ['chr1', 'chr3', 'chr5', 'chr7', 'chr9']
CV_FOLDS = 5
RANDOM_STATE = 42

"""
Functions:
- make_col_combinations()
- train_main()
- compute_auprc()
- plot_precision_recall_curve()
- evaluate_and_save_pr_curves()
"""

def make_col_combinations(df: pd.DataFrame) -> dict[str, list[str]]:
    """
    Generate feature combinations based on the DataFrame columns.
    """
    ensemble_scores = [
        (name, [c for c in df.columns if c.startswith(prefix)]) for name, prefix in ENSEMBLE_SCORE_COLS.items()
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
        if name in RAW_SCORE_COLS:
            combinations[name + "_raw"] = RAW_SCORE_COLS[name]

    return combinations

def train_main(save_dir: Path, df: pd.DataFrame, test_chroms: list[str], log_dir: Path, run_name: str):
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
        e) Compare with using the "raw" feature scores (e.g. SpliceAI_DS) as the only feature,
            to see if the model is actually learning something useful.
    """
    save_dir.mkdir(parents=True, exist_ok=True)
    log_dir.mkdir(parents=True, exist_ok=True)

    df = preprocess(train_mode=True, df=df)
    combinations = make_col_combinations(df)

    # Global train/val/test split before filtering out unscorable rows. This avoids data leakage
    # and allows fair comparison of feature sets with different numbers of unscorable rows.
    # Unscorable rows will be filtered out separately for each feature set during training.
    if TARGET_COL not in df.columns:
        raise ValueError(f"Target column '{TARGET_COL}' is missing after preprocessing. Available columns: {df.columns}")
    if CHROM_COL not in df.columns:
        raise ValueError(f"Chromosome column '{CHROM_COL}' is missing after preprocessing. Available columns: {df.columns}")
    y_all = pd.to_numeric(df[TARGET_COL], errors='raise').astype(int)
    if not set(y_all.unique()).issubset({0, 1}):
        raise ValueError(f"Target column '{TARGET_COL}' must be binary (0 and 1). Found values: {y_all.unique()}")

    X_all = df.drop(columns=[TARGET_COL])
    test_mask = df[CHROM_COL].isin(test_chroms)
    if test_mask.sum() == 0:
        raise ValueError(f"No rows found for test chromosomes {test_chroms}. Check if chromosome naming (e.g. chr1 vs 1).")
    if (~test_mask).sum() == 0:
        raise ValueError(f"No rows remain for train/validation after applying test chromosome split.")
    
    X_train_full = X_all[~test_mask]
    y_train_full = y_all[~test_mask]
    groups_train = df[CHROM_COL].loc[~test_mask]

    X_test_full = X_all.loc[test_mask]
    y_test_full = y_all.loc[test_mask]
    groups_test = df[CHROM_COL].loc[test_mask]

    overlap_chroms = set(groups_train.unique()) & set(groups_test.unique())
    if len(overlap_chroms) != 0:
        raise ValueError(f"Data leakage risk: Chromosomes {', '.join(overlap_chroms)} appear in both train and test sets.")
    
    if groups_train.nunique() < CV_FOLDS:
        raise ValueError(f"Need at least {CV_FOLDS} non-test chromosomes for GroupKFold,"
                         f" but only found {groups_train.nunique()}.")
    
    cv = GroupKFold(n_splits=CV_FOLDS)
    cv_splits = list(cv.split(X_train_full, y_train_full, groups=groups_train))

    results: dict[str, dict] = {} # run name -> stats dict (train auprc, val auprc, etc.)
    run_register: dict[str, dict] = {} # run name -> dict of info to save about that run (e.g. model params, feature columns used, etc.)
    skipped_runs: dict[str, str] = {}  # run name -> reason for skipping
    # # Store models and test sets to evaluate the best one later
    # saved_models = {}
    # saved_test_data = {}
    for feature_group, cols in combinations.items():
        missing_cols = [col for col in cols if col not in X_all.columns]
        error = None
        if missing_cols:
            reason = f"The following columns for feature group '{feature_group}' are missing: {missing_cols}"
            print(reason)
            skipped_runs[feature_group] = reason
            continue

        print(f"Training with feature set: {feature_group} ({len(cols)} features)")

        models = [
            ("RawScoreOnly", (None, {}))
        ] if feature_group.endswith("_raw") else [
            (name, (cls, kwargs)) for name, (cls, kwargs) in CLASSIFIERS.items() if cls is not None
        ]

        for model_name, model_payload in models:
            ModelClass, model_params = model_payload
            run_key = f"{model_name}_{feature_group}"
            fold_train_auprcs: list[float] = []
            fold_val_auprcs: list[float] = []

            for fold_idx, (train_idx, val_idx) in enumerate(cv_splits):
                X_train_fold, y_train_fold = X_train_full.iloc[train_idx], y_train_full.iloc[train_idx]
                X_val_fold, y_val_fold = X_train_full.iloc[val_idx], y_train_full.iloc[val_idx]

                X_train, y_train = filter_bad_features(X_train_fold, y_train_fold, cols)
                X_val, y_val = filter_bad_features(X_val_fold, y_val_fold, cols)

                if len(X_train) == 0 or len(X_val) == 0:
                    error = f"Fold {fold_idx}: No valid training samples after filtering for feature group '{feature_group}'"
                    break

                if ModelClass is None:
                    # Do not fit a model; the prediction is just the max of all "abs(raw score)" from specified columns
                    y_train_preds = X_train.abs().max(axis=1)
                    y_val_preds = X_val.abs().max(axis=1)
                else:
                    clf = ModelClass(random_state=RANDOM_STATE, **model_params)
                    clf.fit(X_train, y_train)

                    y_train_preds = clf.predict_proba(X_train)[:, 1]
                    y_val_preds = clf.predict_proba(X_val)[:, 1]

                fold_train_auprcs.append(compute_auprc(y_train, y_train_preds))
                fold_val_auprcs.append(compute_auprc(y_val, y_val_preds))

            if error is not None:
                print(f"Skipping model {run_key} due to error: {error}")
                skipped_runs[run_key] = error
                continue

            results[run_key] = {
                "train auprc": round(100 * float(np.mean(fold_train_auprcs)), 4),
                "val auprc":   round(100 * float(np.mean(fold_val_auprcs)), 4),
                "train auprc std": round(100 * float(np.std(fold_train_auprcs)), 4),
                "val auprc std":   round(100 * float(np.std(fold_val_auprcs)), 4),
                "num folds": CV_FOLDS,
                "all fold train auprcs": [round(100 * float(auprc), 4) for auprc in fold_train_auprcs],
                "all fold val auprcs": [round(100 * float(auprc), 4) for auprc in fold_val_auprcs],
            }

            run_register[run_key] = {
                "model_name": model_name,
                "feature_group": feature_group,
                "feature_columns": cols,
                "is_raw": feature_group.endswith("_raw"),
                "model_class": ModelClass,
                "model_params": model_params
            }

        # if feature_group.endswith("_raw"):
        #     model_name = "RawScoreOnly"

        #     # Do not fit a model; the prediction is just the max of all "abs(raw score)" from specified columns
        #     y_train_preds = X_train.abs().max(axis=1)
        #     y_val_preds = X_val.abs().max(axis=1)

        #     evaluate_and_save_pr_curves(
        #         y_train, y_train_preds, y_val, y_val_preds,
        #         model_name, feature_group, log_dir, results
        #     )

        #     # Save state for best model evaluation later
        #     run_key = f"{model_name}_{feature_group}"
        #     saved_models[run_key] = None  # No model, just the raw scores
        #     saved_test_data[run_key] = (X_test, y_test)
        #     continue

        # for model_name, (ModelClass, param_grid) in CLASSIFIERS.items():
        #     # Just init with default params for now, can add hyperparameter tuning later
        #     clf = ModelClass()
        #     clf.fit(X_train, y_train)

        #     y_train_preds = clf.predict_proba(X_train)[:, 1]
        #     y_val_preds = clf.predict_proba(X_val)[:, 1]

        #     evaluate_and_save_pr_curves(
        #         y_train, y_train_preds, y_val, y_val_preds,
        #         model_name, feature_group, log_dir, results
        #     )

        #     # Save state for test set evaluation
        #     run_key = f"{model_name}_{feature_group}"
        #     saved_models[run_key] = clf
        #     saved_test_data[run_key] = (X_test, y_test)

    if not results:
        raise ValueError(f"No models were successfully trained. Skipped runs: {skipped_runs}")

    # get best validation AUPRC model
    model_run_keys = [key for key in results.keys() if not key.endswith("_raw")]
    best_model_key = max(model_run_keys, key=lambda run_key: results[run_key]['val auprc'])
    print(f"Best model: {best_model_key} with train AUPRC={results[best_model_key]['train auprc']:.4f} and val AUPRC={results[best_model_key]['val auprc']:.4f}")

    raw_run_keys = [key for key in results.keys() if key.endswith("_raw")]
    best_raw_key = max(raw_run_keys, key=lambda run_key: results[run_key]['val auprc'])
    print(f"Best raw score only model: {best_raw_key} with val AUPRC={results[best_raw_key]['val auprc']:.4f}")

    def evalute_on_test(run_key: str):
        meta = run_register[run_key]
        feature_cols = meta['feature_columns']
        feature_group = meta['feature_group']
        model_name = str(meta['model_name'])

        X_train, y_train = filter_bad_features(X_train_full, y_train_full, feature_cols)
        X_test, y_test = filter_bad_features(X_test_full, y_test_full, feature_cols)

        if len(X_train) == 0 or len(X_test) == 0:
            raise ValueError(f"Skipping test evaluation for {run_key} due to no valid samples after filtering.")

        ModelClass = meta['model_class']
        if ModelClass is None:
            trained_model = None
            y_test_preds = X_test.abs().max(axis=1)
        else:
            model_params = meta['model_params']
            trained_model = ModelClass(random_state=RANDOM_STATE, **model_params)
            trained_model.fit(X_train, y_train)
            y_test_preds = trained_model.predict_proba(X_test)[:, 1]
        
        test_fig, test_auprc = plot_precision_recall_curve(y_test, y_test_preds, model_name, feature_group)
        test_fig.savefig(log_dir / f"{run_key}_{run_name}_test_pr_curve.png")
        plt.close(test_fig)
        return trained_model, round(100 * test_auprc, 4)

    # evaluate the best model on the test set and save metrics and PR curve
    best_model, best_model_test_auprc = evalute_on_test(best_model_key)
    assert best_model is not None, "Best model should not be None since it was trained above. Check for errors during training."
    results[best_model_key]["test auprc"] = best_model_test_auprc
    print(f"Best model test AUPRC: {best_model_test_auprc:.4f}")

    _, best_raw_test_auprc = evalute_on_test(best_raw_key)
    results[best_raw_key]["test auprc"] = best_raw_test_auprc
    print(f"Best raw-score-only model test AUPRC: {best_raw_test_auprc:.4f}")

    # Save all scores, and the best model and according columns
    with open(save_dir / f"{best_model_key}_{run_name}_columns.json", 'w') as f:
        json.dump(run_register[best_model_key]['feature_columns'], f, indent=4)
    with open(save_dir / f"{best_model_key}_{run_name}_model.pkl", 'wb') as f:
        pickle.dump(best_model, f)

    best_model_summary = {
        "description": "best model by validation AUPRC among all non-raw-score-only models",
        "model_key": best_model_key,
        "validation_auprc": results[best_model_key]['val auprc'],
        "test_auprc": results[best_model_key]['test auprc'],
        "params": run_register[best_model_key]['model_params']
    }

    raw_summary = {
        "description": "best raw-score-only model by validation AUPRC",
        "model_key": best_raw_key,
        "validation_auprc": results[best_raw_key]['val auprc'],
        "test_auprc": results[best_raw_key]['test auprc'],
    }

    # Save as list sorted by desc order of val AUPRC
    sorted_results = dict(sorted(results.items(), key=lambda item: item[1]['val auprc'], reverse=True))
    payload = {
        "best_model": best_model_summary,
        "best_raw_score_model": raw_summary,
        "all_results": sorted_results,
        "skipped_runs": skipped_runs,
    }

    with open(log_dir / f"training_results_{run_name}.json", "w") as f:
        json.dump(payload, f, indent=4)

def compute_auprc(y_true, y_scores) -> float:
    precision, recall, _ = precision_recall_curve(y_true, y_scores)
    return float(auc(recall, precision))

def plot_precision_recall_curve(y_true, y_scores, model_name, feature_set_name) -> tuple[Figure, float]:
    """
    Returns figure object.
    """
    precision, recall, thresholds = precision_recall_curve(y_true, y_scores)

    # F1 for all threshold points (same length as thresholds)
    # precision_recall_curve returns one extra point where recall=0;
    # thresholds has length len(precision) - 1, so slice.
    f1_scores = 2 * precision[:-1] * recall[:-1] / (precision[:-1] + recall[:-1] + 1e-12)

    best_idx = np.argmax(f1_scores)
    best_f1 = f1_scores[best_idx]
    best_threshold = thresholds[best_idx]

    pr_auc = float(auc(recall, precision))

    # Use Object-Oriented API
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(recall, precision, label=f'PR curve (AUC={pr_auc:.3f})')

    # Mark the best F1 point
    ax.scatter(recall[best_idx], precision[best_idx], color='red', label=f'Best F1={best_f1:.3f}')

    ax.set_xlabel('Recall')
    ax.set_ylabel('Precision')
    ax.set_title(f'PR Curve: {model_name}\nFeatures: {feature_set_name}')
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.05) # Small padding at top
    ax.legend(loc='lower left')
    ax.grid(True, linestyle='--', alpha=0.7)

    return fig, pr_auc

def evaluate_and_save_pr_curves(y_train, y_train_preds, y_val, y_val_preds,
                                model_name: str, feature_group: str, log_dir: Path):
    """Helper function to plot, save, and record Precision-Recall curves."""
    train_fig, train_auprc = plot_precision_recall_curve(y_train, y_train_preds, model_name, feature_group)
    val_fig, val_auprc = plot_precision_recall_curve(y_val, y_val_preds, model_name, feature_group)

    # results[f"{model_name}_{feature_group}"] = {
    #     "train auprc": round(100 * train_auprc, 4),
    #     "val auprc": round(100 * val_auprc, 4)
    # }

    # Save the figures
    train_fig.savefig(log_dir / f"{model_name}_{feature_group}_train_pr_curve.png")
    val_fig.savefig(log_dir / f"{model_name}_{feature_group}_val_pr_curve.png")

    # Close figures to free up memory
    plt.close(train_fig)
    plt.close(val_fig)

