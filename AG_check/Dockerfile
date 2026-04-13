FROM python:3.10-slim
WORKDIR /app

# Update package lists and install necessary dependencies
RUN apt-get update && \
    apt-get install -y \
        autoconf \
        automake \
        gcc \
        make \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        curl \
        bzip2 \
        git \
        procps && \
    rm -rf /var/lib/apt/lists/*

# Download and install HTSlib (which contains bgzip)
RUN curl -L https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 -o htslib.tar.bz2 && \
    tar xvjf htslib.tar.bz2 && \
    cd htslib-1.17 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf htslib-1.17 htslib.tar.bz2

# Download and install bcftools
RUN curl -L https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 -o bcftools.tar.bz2 && \
    tar xvjf bcftools.tar.bz2 && \
    cd bcftools-1.17 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -rf bcftools-1.17 bcftools.tar.bz2

RUN pip install --upgrade pip --no-cache-dir && \
    pip install --no-cache-dir pysam Biopython pandas vcf2pandas scikit-learn matplotlib

RUN curl -L https://github.com/brentp/vcfanno/releases/download/v0.3.5/vcfanno_linux64 -o /usr/local/bin/vcfanno && \
    chmod +x /usr/local/bin/vcfanno

LABEL version='v1.6'
