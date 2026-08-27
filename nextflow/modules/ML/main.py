import argparse
from pathlib import Path

from vcf2pandas import vcf2pandas

from pipeline_constants import (
    INFO_FIELDS,
)

from train import train_main
from infer import infer_main

def existing_file(path_str: str) -> Path:
    path = Path(path_str)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"Expected an existing file path, got: {path_str}")
    return path

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Train or run inference for the Introme ensemble score model")
    subparsers = parser.add_subparsers(dest="mode", required=True, help="Mode of operation: 'train' or 'infer'")

    infer_parser = subparsers.add_parser("infer", help="Run inference using a saved model and feature columns")
    infer_parser.add_argument("--model-path", required=True, type=existing_file, help="Path to an existing model pickle (.pkl) file")
    infer_parser.add_argument("--columns-path", required=True, type=existing_file, help="Path to an existing columns JSON file")
    infer_parser.add_argument("--input-vcf", required=True, type=existing_file, help="Path to the input VCF file with splicing annotations")
    infer_parser.add_argument("--output-tsv", required=True, type=Path, help="Path to write the output TSV scores")

    train_parser = subparsers.add_parser("train", help="Train candidate models and save selected model + columns")
    train_parser.add_argument("--save-dir", required=True, type=Path, help="Directory to save the best model pickle and columns JSON")
    train_parser.add_argument("--input-vcf", required=True, type=existing_file, help="Path to the input VCF file with splicing annotations and labels")
    train_parser.add_argument("--log-dir", required=True, type=Path, help="Directory to save training logs, metrics, and plots")
    train_parser.add_argument("--run-name", required=True, type=str, help="Name for this training run (used in saved files)")
    train_parser.add_argument("--test-chroms", required=True, nargs='+', help=f"List of chromosomes to use as test set")

    return parser

if __name__ == "__main__":
    parser = build_parser()
    args = parser.parse_args()

    df = vcf2pandas(str(args.input_vcf),
        remove_empty_columns=False,
        info_fields=INFO_FIELDS
    )

    if (args.mode == "train"):
        train_main(
            save_dir=args.save_dir,
            df=df,
            test_chroms=args.test_chroms,
            log_dir=args.log_dir,
            run_name=args.run_name
        )
    else:
        infer_main(
            model_path=args.model_path,
            df=df,
            columns_path=args.columns_path,
            output_path=args.output_tsv
        )

