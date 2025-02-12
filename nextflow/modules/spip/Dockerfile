# Use the r-base image on amd64 architecture
FROM --platform=linux/amd64 r-base:latest

# Suppress interactive prompts from apt
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies needed for SPiP (git, wget, etc.)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        wget \
    && rm -rf /var/lib/apt/lists/*

# Install the required R packages
RUN R -e "install.packages(c('foreach', 'doParallel', 'randomForest'), \
          repos='http://cran.us.r-project.org')"

# Define build-time variables
ARG SPIP_REPO=https://github.com/LBGC-CFB/SPiP.git
ARG SPIP_VERSION=cae95fe0ee7a2602630b7a4eacbf7cfa0e1763f0

# Clone the SPiP repository and download reference files
RUN git init SPiP \
    && cd SPiP \
    && git remote add origin ${SPIP_REPO} \
    && git fetch --depth 1 origin ${SPIP_VERSION} \
    && git checkout FETCH_HEAD \
    && cd ./RefFiles \
    && wget -O transcriptome_hg38.RData \
      https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg38.RData/download \
    && wget -O transcriptome_hg19.RData \
      https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg19.RData/download

# Set SPiP as the working directory
WORKDIR /SPiP

LABEL version="v1.1"
