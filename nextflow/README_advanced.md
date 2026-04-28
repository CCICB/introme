# Introme v2: Advanced Guide

This guide covers architecture, module-level design, assets, containers, model training, and migration details.

## Table of Contents

- Basic setup and run instructions: [README.md](README.md)
- This file contains deep technical documentation for development and extension.

1. [Module Map](#module-map)
1. [Legacy v1 -> v2 Step Mapping](#legacy-v1---v2-step-mapping)
1. [Annotation and Asset Files](#annotation-and-asset-files)
1. [Container and Docker Layout](#container-and-docker-layout)
1. [Extending Introme](#extending-introme)
1. [How to Train the Consensus Scoring Model](#how-to-train-the-consensus-scoring-model)
1. [Known Caveats and TODOs](#known-caveats-and-todos)
1. [Recommended Operational Workflow](#recommended-operational-workflow)

## Module Map

Each module script is in `modules/MODULE.nf` or `modules/MODULE/MODULE.nf` for containerised tools.
Modules ingest defined inputs, and upon completion, respectively expose selected files via symlink in their designated `work/MODULE/` folder.

| Stage | Module | Key inputs | Assets/static dependencies | Key outputs |
| --- | --- | --- | --- | --- |
| Subset + normalise | `data_preprocessing` | Input VCF, reference FASTA, input GTF | `chrRename.tsv` | `*.subset.vcf.gz`, `sorted.gtf.gz` |
| Quality filter | `quality_filter` | Input VCF, params (`quality_filter`, `min_QUAL`, `min_DP`, `min_AD`) | - | `*.quality_filter.vcf.gz` |
| Pre-annotation | `variant_info` | Filtered/unfilterd VCF, sorted GTF | vcfanno `assets/` files: `conf_pre_anno.lua`, `tomls/gencode.[build].toml` | `*.variant_info.vcf.gz` and stripped + rmanno versions |
| | | | | |
| SpliceAI | `spliceai` | rmanno VCF, FASTA | SpliceAI anno DB (`params.spliceai_db`, downloaded at runtime) | `*.spliceai.vcf` |
| MMSplice | `mmsplice` | rmanno VCF, FASTA, GTF | `modules/mmsplice/run_mmsplice.py` | `*.mmsplice.vcf` |
| Pangolin | `pangolin` | rmanno VCF, FASTA | Pangolin anno DB (`params.pangolin_db`, downloaded at runtime) | `*.pangolin.vcf` |
| SPiP | `spip` | rmanno VCF | SPiP resources are included in container image | `*.spip.vcf` |
| Spliceogen | `spliceogen` | rmanno VCF, FASTA, GTF | - | `*.spliceogen.tsv` |
| Introme feature modules | `introme_functions` | stripped VCF, FASTA | `AG_check/AG_check.py`, `ESE/scoring.py`, `assets/introme_annotate.vcf` | `*.ag_check.vcf.gz`, `*.ESE.tsv.gz` |
| | | | | |
| Ensemble collation | `splicing_anno` | `variant_info` VCF, tool score and `.vcf`/`.tsv` outputs | vcfanno `assets/` files: `conf_ensemble.lua`, `tomls/annotate.[build].toml`, `tomls/vcfanno_splicing_ensemble.toml`, `branchpointer/[build].bed.gz`, `regions/[build].bed.gz`, `U12/[build].bed.gz` | `*.ensemblescored.vcf.gz` |
| ML inference | `ensemble_infer` | Scored VCF from `splicing_anno` | `ESE/ML/main.py`, model `.pkl`, columns `.json` | `*.introme.predictions.tsv` |
| ML training | `ensemble_train` | Scored VCF from `splicing_anno` | `ESE/ML/main.py`, test chromosome set, train config params (`params.json`) | `models/*`, `logs/*` |

## Legacy v1 -> v2 Step Mapping

| v1 (`run_introme.sh`) | v2 module/process |
| --- | --- |
| Step 1: subset VCF | `data_preprocessing` |
| Step 2: quality filter | `quality_filter` (optional) |
| Step 3: annotate VCF | `variant_info` |
| Step 4: annotation-value filtering based on gnomAD_PM_AF | Disabled: Commented out AF anno in .toml, and vcfanno in `variant_info` |
| Step 5: MMSplice/SpliceAI (+tool stage) | `spliceai`, `mmsplice`, `pangolin`, `spip`, `spliceogen` |
| Step 6: Introme functions | `introme_functions` + `splicing_anno` |
| Step 7: TSV export | produced via ensemble inference output TSV |
| Step 8: consensus scoring | `ensemble_infer` or `ensemble_train` |

## Annotation and Asset Files

### Lua

- `assets/conf_pre_anno.lua`: helper functions for pre-annotation (gene, strand, region logic).
- `assets/conf_ensemble.lua`: parsing/conversion of output formats from SpliceAI/MMSplice/Pangolin/SPiP for vcfanno.

### TOML

- `assets/tomls/gencode.hg[19|38].toml`: pre-annotation sources for gene-level context.
- `assets/tomls/annotate.hg[19|38].toml`: branchpointer/regions/U12 annotations.
- `assets/tomls/vcfanno_splicing_ensemble.toml`: tool-score ingestion and INFO mapping for the final scored VCF.

### Build-specific annotation datasets

- `assets/branchpointer/*.bed.gz(.tbi)`
- `assets/regions/*.bed.gz(.tbi)`
- `assets/U12/*.bed.gz(.tbi)`

### Utility/reference assets

- `assets/chrRename.tsv`: chromosome renaming map used early in preprocessing.
- `assets/introme_annotate.vcf`: template/header utility for AG/ESE related outputs.
- `assets/models/`: default model artifacts used by inference mode. Contains `.pkl` model and corresponding `.json` columns file.

## Container and Docker Layout

Runtime images are configured via `params.json` and selected by profile in `nextflow.config`.

Tool image params include:

- Images with cpu/gpu verisons
  - `spliceai_[cpu|gpu]_docker_container`
  - `mmsplice_[cpu|gpu]_docker_container`
  - `pangolin_[cpu|gpu]_docker_container`
- `spip_docker_container`, `spliceogen_docker_container`
- `data_preprocessing_docker_container`, `variant_info_docker_container`, `introme_functions_docker_container`

Source Dockerfiles for maintenance and rebuilds:

- `modules/spliceai/Dockerfile.[cpu|gpu]`
- `modules/mmsplice/Dockerfile.[cpu|gpu]`
- `modules/pangolin/Dockerfile.[cpu|gpu]`
- `modules/spip/Dockerfile`
- `modules/spliceogen/Dockerfile`
- `../AG_check/Dockerfile`
- Dockerfile sources for data_preprocessing, variant_info are MIA but can be reverse enginerred by inspecting "Layers" of image.

> NOTE: GPU images are large (each 10GB) due to needing GPU libraries. SpliceAI and MMSplice use tensorflow. Pangolin uses Pytorch. Images can probably be slimmed down, probably also can merge the tensorflow images.

## Extending Introme

This section describes how to introduce or remove features in two common scenarios.

1. Case A: a new external tool is called and dynamically produces scores.
1. Case B: an additional annotation source is added via TOML only.

### Case A: Tool-calling integration (new score-producing tool)

#### Step 1. Dockerfile and process wiring

1. Create or update container build files under `modules/TOOL_NAME/`, including `Dockerfile` or `Dockerfile.cpu` and `Dockerfile.gpu` plus optional helper runner scripts.
1. Add a process module `modules/TOOL_NAME/TOOL_NAME.nf` that consumes the expected VCF/FASTA/GTF inputs, emits an output file (VCF/TSV), and publishes outputs to `output/TOOL_NAME`.
1. Wire the process into orchestration in `main.nf`: include the module, call it in workflow order (usually after `variant_info`), and pass its output into `splicing_anno`.
1. Register runtime image/config knobs in `params.json` (container names, db filenames, toggles) and `nextflow.config` profiles (cpu/mem/gpu label if needed).

#### Step 2. Parse tool output into ensemble features

1. Add tool extraction logic to `assets/tomls/vcfanno_splicing_ensemble.toml` with `[[annotation]]` for raw field ingestion, `[[postannotation]]` for parsing and feature extraction, and `op = "delete"` cleanup for temporary fields.
1. If the tool output format is complex, add a parser function in `assets/conf_ensemble.lua` and call it from TOML.
1. Ensure resulting INFO names are stable and model-friendly (typically prefixed by tool name).
1. Update all relevant constants in `ESE/ML/pipeline_constants.py`:
    - `INFO_FIELDS` (maps incoming INFO fields into ML naming),
    - `ENSEMBLE_SCORE_COLS` (group definition used for feature combinations),
    - `RAW_SCORE_COLS` (raw-score-only benchmark comparisons).
    - When adding/removing categorical one-hot fields, register it to `DUMMY_COLS`, list out all possible values cf. `GENE_REGIONS`/`VARIANT_TYPE`, and add it to `columns` keyword of `pd.get_dummies()` in `preprocess()`.

> Note: it is possible to add a feature but not train/infer on it, i.e. just to report it. Simply do not add it to `ENSEMBLE_SCORE_COLS`/`RAW_SCORE_COLS`.

To remove a tool, remove its blocks from `vcfanno_splicing_ensemble.toml` and remove its entries from the Python constants.

#### Step 3. Missing-feature imputation and tuple sanitisation

1. Update preprocessing in `ESE/ML/utils.py` inside `preprocess(...)`, where `assert_and_convert_single_float_tuples_allow_dot(...)` is called.
1. Add or remove tool prefixes in the imputation controls: `dots_largeminus`, `dots_zeros`, `nans_zeros`, and `nans_large_minus`.
1. Choose strategy intentionally: use large negative sentinels for unscorable or missing-is-informative values, and use zeros for neutral absence.
1. Verify downstream filtering behavior in `filter_bad_features(...)` remains aligned with your missingness policy.

### Case B: Annotation-only integration (no new tool call)

This path is for adding static annotations (for example branchpointer-like tracks or extra BED/VCF-based fields).

#### Step 1. Annotation TOML changes

1. Add annotation entries in `assets/tomls/annotate.hg38.toml` (and usually `annotate.hg19.toml` for parity), including `[[annotation]]` sources, optional `[[postannotation]]` transforms, and optional cleanup with `delete`.
1. Ensure files are present under `assets/` and referenced with the same relative paths expected by `vcfanno`.
1. No new Dockerfile is required if this is purely annotation-file driven.

Step 2-3: Follow the above method for Case A.

### Step 4. Training and evaluation (applies to both cases)

1. Run training through Nextflow (`--ml_mode train`) after every feature-schema change.
1. Compare best model versus raw-score-only baselines in training outputs (`training_results_*.json`, PR curves).
1. Confirm exported model/columns artifacts are regenerated together and versioned as a pair.
1. Re-run inference with the new model and validate TSV output columns against expectations.
1. If removing tools/features, retrain before production use to avoid stale expected columns.

## How to Train the Consensus Scoring Model

Training uses `../ESE/ML/main.py train` through `ensemble_train`.

Required training params:

1. `--ml_mode train`
1. `--ml_save_dir`
1. `--ml_log_dir`
1. `--ml_run_name`

Optional params:

1. `--ml_test_chroms` (defaults to: `chr1 chr3 chr5 chr7 chr9`)

Example command:

```bash
nextflow run . \
  --vcf /path/to/input.vcf.gz \
  --ref_genome /path/to/hg38.fa \
  --gtf /path/to/gencode.v46.annotation.gtf.gz \
  --chrRename ./assets/chrRename.tsv \
  --prefix Apr2026 \
  --quality_filter false \
  --ml_mode train \
  --ml_save_dir ./output/ensemble_train/models \
  --ml_log_dir ./output/ensemble_train/logs \
  --ml_run_name mytrainingrun \
  --ml_test_chroms chr1,chr3,chr5,chr7,chr9 \
  -params-file ./params.json \
  -profile gpu \
  -resume
```

Training artifacts are published to:

- models: `${ml_save_dir}/models/*`
- logs: `${ml_log_dir}/logs/*`

## Known Caveats and TODOs

### Relating to Splice Prediction Tools

1. Some variants are skipped by certain tools, resulting in the need for imputing/flagging missingness for the ensemble model. This is transparent; variants are kept through to report generation (although Introme scores are currently appended to the imputed dataframe, not the version with missing inputs). Some fail cases I've observed:
    - Spliceogen in all but 1 in 25k variants has a missing feature
    - Variant not in gene region (not sure whether it's according to tool's own gene annotation, or provided one)
    - Variant too close to chrom ends
    - Pangolin doesn't support `(len(ref) != 1 and len(alt) != 1 and len(ref) != len(alt))` ([Github source](https://github.com/tkzeng/Pangolin/blob/5cf94b8db938c658391b4305cd7ce33297d44ff7/pangolin/pangolin.py#L88))
    - Lowercase variant bases, e.g. `chr11:47334593 ctCTG>c`
    - For further investigation, missing features can be filtered for by feature name in the jupyter notebook.
1. SpliceAI, Pangolin script entrypoints are fine for individual variant calls, but especially for GPU-enabled machines __do not take advantage of batched processing__. It should not take 90 minutes to process 25k variants. For comparison, training OpenSpliceAI from scratch takes 2 hours on the same hardware.
1. The Pangolin module will download its 300MB annotation file each time it runs -- SpliceAI's 10MB is not that bad. Since this is not removed, these files will accumulate on your computer in the work directory. Even if it were to be removed, isn't it wasteful to host a server to download this to you? TODO: either have an ability to specify a Pangolin annotation file path as main params, or dynamically create the annotation database each time, as per their instructions:

    - > Create an annotation database from a GTF file using scripts/create_db.py. This will take several minutes. ([Github source](https://github.com/tkzeng/Pangolin/blob/main/README.md))

1. SQUIRL integration exists (`modules/squirl.nf`) but is currently not active in `main.nf`. Some other Dockerfile sources e.g. Absplice are being kept in Confluence.
1. `MNV.sh` scoring logic is deprecated, as we rely on the implementation in bw2 fork, introduced in Apr 2023 ([Github source](https://github.com/bw2/SpliceAI/commit/1dcd441d4e931909007f06a60c9e285c994699b0)).

### Relating to Introme ML/functions

1. Hyperparameter tuning for the current training path remains a TODO. Current behaviour is defaulting to sklearn’s default initialisation params. For instance, this creates quite deep trees in RF (which increases the model size, probably prone to overfit too).
1. Train/inference step in nextflow is set to "no cache" because otherwise -- outside of changes to the main.py script -- if you change the script, it would not detect changes and will happily reuse past results. Another solution is to import the path of the entire ml source directory.
1. Dual +/- gene strand is not supported in ESE/scoring.py (scoring goes to + only case). Double check if this is indeed the case. +/- cases should be split into two records for simplicity.

### Others

- To reiterate, hg19 has not been tested, and doubtful to work until annotation `.toml`'s are updated.
- Current `main.nf` quality filter call path bug was corrected, but not tested.
- Using nf-schema to validate params
- Tips for MLOps? Experiment tracking, model versioning etc.
- CI/CD for automated tests?

## Recommended Operational Workflow

1. Validate asset TOML paths and build-specific annotations (`hg19` vs `hg38`) before first run.
1. Run a small test VCF in inference mode with `-resume` to validate environment and containers.
1. Scale to full cohort runs.
1. Use training mode to refresh model artifacts when updating features/tool versions.
