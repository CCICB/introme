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

# Clone the SPiP repository and download reference files
RUN git clone https://github.com/raphaelleman/SPiP && \
    cd SPiP/RefFiles && \
    wget -O transcriptome_hg38.RData \
      https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg38.RData/download && \
    wget -O transcriptome_hg19.RData \
      https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg19.RData/download

# Set SPiP as the working directory
WORKDIR /SPiP
