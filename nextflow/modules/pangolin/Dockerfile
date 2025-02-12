# Use the Miniconda base image
FROM continuumio/miniconda3:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive \
LANG=C.UTF-8 \
LC_ALL=C.UTF-8 \
PATH=/opt/conda/bin:$PATH

# Install system packages and clean up
RUN apt-get update && \
apt-get install -y --no-install-recommends apt-utils wget git && \
apt-get -y upgrade && \
rm -rf /var/lib/apt/lists/*

# Update Conda and install Python and packages
RUN conda update -n base -c defaults conda && \
conda install -y python=3.9 && \
conda config --add channels bioconda && \
conda install -y -c conda-forge pyvcf && \
conda install -y pytorch cpuonly -c pytorch && \
conda clean -afy

# Define build-time variable for Pangolin commit
ARG PANGOLIN_COMMIT=fd12887eff04720ea0b4cdd02373f23132b39774

# Upgrade pip and install Python packages
RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir gffutils biopython pandas pyfastx && \
    pip install --no-cache-dir pangolin@git+https://github.com/bw2/Pangolin.git@${PANGOLIN_COMMIT}

LABEL version="v1.2"
