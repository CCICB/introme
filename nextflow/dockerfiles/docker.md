# Docker Image Changelog and Build Commands

## Changelog

### MMSplice

3.0.1: GPU image builds torch CUDA dependencies through pip, not as base image.
3.0.0: Using anaconda/miniconda:26.9.4 as base for both CPU. Rich labels.

### Pangolin

2.0.1: procps back in the base images. GPU image builds torch CUDA dependencies through pip, not as base image.
2.0.0: Using anaconda/miniconda:26.9.4 as base for CPU. Rich labels. procps tacked on in annotation image in version b. Annotation image bakes in pangolin generated db annotations from provided annotation gtf.

### SpliceAI

2.0.0: Using anaconda/miniconda:26.9.4 as base for CPU. Rich labels. Annotation image bakes in bw-spliceai generated db annotations from provided annotation gtf.

### utils

1.0.1: added validation of pinned versions and image labels.
1.0.0: introme-utils consolidates previous data_preprocessing and variant_info dockers into one versions are as follows:

- bcftools 1.13
- SAMtools <unused; removed>
- BEDTools v2.30.0
- Tabix 1.13+ds
- BGZip 1.13+ds (same version as tabix)
- vcfanno version 0.3.3 [built with go1.18.1]

## Build Commands

```sh
# mmsplice
# cd into nextflow/dockerfiles/mmsplice
# build script refers to run_mmsplice.py
docker build -f Dockerfile.cpu -t introme-mmsplice:3.0.0-runtime-cpu .
docker build -f Dockerfile.gpu -t introme-mmsplice:3.0.0-runtime-gpu .

# Pangolin
# Two stage build: 1. cpu/gpu base 2. bake in pangolin generated db from provided annotation gtf
docker build -f nextflow/dockerfiles/pangolin/Dockerfile.cpu  \
    --build-arg IMAGE_VERSION=2.0.0 \
    -t introme-pangolin:2.0.0-runtime-cpu     .

docker build -f nextflow/dockerfiles/pangolin/Dockerfile.gpu \
    --build-arg IMAGE_VERSION=2.0.0 \
    -t introme-pangolin:2.0.0-runtime-gpu     .

docker build -f nextflow/dockerfiles/pangolin/Dockerfile.annotation \
    --build-arg BASE_IMAGE=introme-pangolin:2.0.0-runtime-cpu \
    --build-arg GENCODE_VERSION=50 \
    --build-arg GENCODE_GTF=</path/to/gtf> \
    -t introme-pangolin:2.0.0-cpu-gencode50     .

docker build -f nextflow/dockerfiles/pangolin/Dockerfile.annotation \
    --build-arg BASE_IMAGE=introme-pangolin:2.0.0-runtime-gpu \
    --build-arg GENCODE_VERSION=50 \
    --build-arg GENCODE_GTF=</path/to/gtf> \
    -t introme-pangolin:2.0.0-gpu-gencode50     .

# spliceai
docker build -f nextflow/dockerfiles/spliceai/Dockerfile.cpu  \
    --build-arg IMAGE_VERSION=2.0.0 \
    -t introme-spliceai:2.0.0-runtime-cpu     .

docker build -f nextflow/dockerfiles/spliceai/Dockerfile.gpu \
    --build-arg IMAGE_VERSION=2.0.0 \
    -t introme-spliceai:2.0.0-runtime-gpu .

docker build -f nextflow/dockerfiles/spliceai/Dockerfile.annotation \
    --build-arg BASE_IMAGE=introme-spliceai:2.0.0-runtime-cpu \
    --build-arg GENCODE_VERSION=50 \
    --build-arg GENCODE_GTF=</path/to/gtf> \
    -t introme-spliceai:2.0.0-cpu-gencode-v50   .

docker build -f nextflow/dockerfiles/spliceai/Dockerfile.annotation \
    --build-arg BASE_IMAGE=introme-spliceai:2.0.0-runtime-gpu \
    --build-arg GENCODE_VERSION=50 \
    --build-arg GENCODE_GTF=</path/to/gtf> \
    -t introme-spliceai:2.0.0-gpu-gencode-v50   .

# utils
docker build -f nextflow/dockerfiles/introme-utils/Dockerfile -t introme-utils:1.0 .

```
