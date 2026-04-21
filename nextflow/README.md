# Introme v2: Basic Guide

Introme predicts the impact of variants on gene splicing by integrating outputs from several splice prediction tools, additional splicing rules, and gene architecture features into a ML model.

**Version 2.0** migrates the original shell workflow from `run_introme.sh` to a modular Nextflow pipeline in `nextflow/main.nf`.

## Table of Contents

- Basic usage and setup are in this file.
- Advanced architecture, module-level design, assets, and model training details are in [README_advanced.md](README_advanced.md).

1. [Changes from v1 to v2](#changes-from-v1-to-v2)
1. [High-Level Architecture](#high-level-architecture)
1. [Setup Guide](#setup-guide)
1. [Configuration Files](#configuration-files)
1. [How to Run Introme 2.0](#how-to-run-introme-20)
1. [Next Steps](#next-steps)

## Changes from v1 to v2

| Area | Introme v1 | Introme v2 (Nextflow) | v2 Advantages |
| --- | --- | --- | --- |
| Architecture | Single shell script (`run_introme.sh`) | Workflow in `main.nf` + modules in `nextflow/modules/` | Easier to maintain and extend |
| Models | Spliceogen, MMSplice, SpliceAI | Spliceogen, MMSplice, SpliceAI, Pangolin, SPiP | Further tools can be incorporated |
| Dependencies | Mixed local installs and container calls | All modules are configured with container images | Simpler environment setup |
| Parallelism | Mostly sequential | Nextflow analyses task dependencies | Automatically managed parallel execution |
| Resume/caching | Manual reruns | `-resume` and work-dir caching | Better fault recovery |
| Hardware profiles | N/A | Customiseable per-module resource profiles | Can be configured to suit availability (e.g. GPU, HPC) |
| Consensus model | R script based scoring (`consensus_scoring.R`) | Python ML entrypoint with `infer` and `train` modes | Unified training/inference path |

## High-Level Architecture

Execution order in `main.nf`:

- `data_preprocessing`: chromosome rename, sort, normalize, subset to GTF or BED regions.
- `quality_filter` (optional): hard quality filtering using `min_QUAL`, `min_DP`, and `min_AD`.
- `variant_info`: pre-annotation with `vcfanno` (`conf_pre_anno.lua` + `gencode.[genome_build].toml`).
- Parallel tool stage: `spliceai`, `mmsplice`, `pangolin`, `spip`, `spliceogen`.
- `introme_functions`: AG/GT gain-loss checks, SNV/INDEL/INSDEL and ESE motif scoring.
- `splicing_anno`: merges all tool outputs and Introme annotations using `vcfanno`.
- Ensemble ML:
  - `ensemble_infer` (default) calculates an Introme score and outputs everything in a `.tsv` report, or
  - `ensemble_train` (`--ml_mode train`): model training on features + evaluation, logging, and selection.

## Setup Guide

### Prerequisites

1. Linux host with enough disk space for images (~12GB for CPU, or ~35GB for GPU), and work and output directories.
1. Nextflow installed and available on `PATH`.
1. Docker installed and running. Docker images will be automatically pulled as needed.
1. For GPU-enabled runs: NVIDIA device and drivers.

### Required input files

1. Input VCF
1. Reference genome FASTA hg19/hg38.
1. Gene annotation GTF.

### Optional input files

1. BED file for region-restricted analysis (`--bed`).
<!-- 1. Custom ML model + columns JSON for inference mode. -->

## Configuration Files

| ! NOTE                      |
|:----------------------------|
| Paths are resolved relative to the terminal (alias `$launchDir`). Use `$projectDir` if referring to a file relative to the folder containing `main.nf`|

- `params.json`:
  - User-editable run inputs, including paths to VCF, FASTA, GTF, and BED files.
  - Values can be overwritten by passing in `--param_name VALUE` in the CLI.
- `update_params.sh`:
  - TODO: fix. Helper script to patch selected `params.json` keys from CLI.
- `nextflow.config`:
  - Contains constants for pipeline behaviour. Add your custom resource profiles here!
  - Profile `standard`: CPU-oriented resource settings.
  - Profile `gpu`: GPU-enabled execution, applies Docker `--gpus all` to all modules with label `gpu` in defined in .nf script.

## How to Run Introme 2.0

It is recommended, but optional, to run from this directory (`introme/nextflow`):

```bash
nextflow run . \
  --vcf /path/to/input.vcf.gz \
  --ref_genome /path/to/hg38.fa \
  --gtf /path/to/gencode.v46.annotation.gtf.gz \
  --prefix my_run \
  --quality_filter false \
  -params-file ./params.json \
  -profile standard \
  -resume
```

For GPU-enabled execution:

```bash
nextflow run . \
  --vcf /path/to/input.vcf.gz \
  --ref_genome /path/to/hg38.fa \
  --gtf /path/to/gencode.v46.annotation.gtf.gz \
  --prefix my_run_gpu \
  --quality_filter false \
  -params-file ./params.json \
  -profile gpu \
  -resume
```

`--param` flags are optional if they are already defined in `params.json`.

Select files from each step are published in `./output` by default, split by process directory.
For detailed output mapping and training artifact locations, see [Output Structure](README_advanced.md#output-structure).

## Next Steps

- If you are developing or extending the pipeline: continue with [README_advanced.md](README_advanced.md).
