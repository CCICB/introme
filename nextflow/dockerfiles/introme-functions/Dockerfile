ARG PYTHON_VERSION=3.10
FROM python:${PYTHON_VERSION}-slim-trixie

ARG IMAGE_VERSION=1.7.0
ARG PYTHON_VERSION

ARG BCFTOOLS_VERSION=1.17
ARG TABIX_VERSION=1.17
ARG VCFANNO_VERSION=0.3.5

# Python package versions
ARG NUMPY_VERSION=2.2.6
ARG BIOPYTHON_VERSION=1.87
ARG PANDAS_VERSION=2.3.3
ARG PYSAM_VERSION=0.22.1
ARG SCIKIT_LEARN_VERSION=1.7.2
ARG SCIPY_VERSION=1.15.3
ARG VCF2PANDAS_VERSION=0.2.0
ARG MATPLOTLIB_VERSION=3.10.8

# Environment configuration
ARG DEBIAN_FRONTEND=noninteractive

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1

WORKDIR /app

# Update package lists and install necessary dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
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

# Install HTSlib (including bgzip and tabix)
RUN curl -fsSL \
        "https://github.com/samtools/htslib/releases/download/${TABIX_VERSION}/htslib-${TABIX_VERSION}.tar.bz2" \
        -o /tmp/htslib.tar.bz2 && \
    tar -xjf /tmp/htslib.tar.bz2 -C /tmp && \
    cd /tmp/htslib-${TABIX_VERSION} && \
    ./configure --prefix=/usr/local && \
    make -j"$(nproc)" && \
    make install && \
    ldconfig && \
    rm -rf /tmp/htslib-${TABIX_VERSION} /tmp/htslib.tar.bz2

# Install BCFtools using the installed HTSlib
RUN curl -fsSL \
        "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" \
        -o /tmp/bcftools.tar.bz2 && \
    tar -xjf /tmp/bcftools.tar.bz2 -C /tmp && \
    cd /tmp/bcftools-${BCFTOOLS_VERSION} && \
    ./configure --prefix=/usr/local --with-htslib=system && \
    make -j"$(nproc)" && \
    make install && \
    rm -rf /tmp/bcftools-${BCFTOOLS_VERSION} /tmp/bcftools.tar.bz2

RUN curl -fsSL \
        "https://github.com/brentp/vcfanno/releases/download/v${VCFANNO_VERSION}/vcfanno_linux64" \
        -o /usr/local/bin/vcfanno && \
    chmod +x /usr/local/bin/vcfanno

# Install Python dependencies
RUN python -m pip install \
        numpy==${NUMPY_VERSION} \
        biopython==${BIOPYTHON_VERSION} \
        pandas==${PANDAS_VERSION} \
        pysam==${PYSAM_VERSION} \
        scikit-learn==${SCIKIT_LEARN_VERSION} \
        scipy==${SCIPY_VERSION} \
        vcf2pandas==${VCF2PANDAS_VERSION} \
        matplotlib==${MATPLOTLIB_VERSION}

# Verify installed bioinformatics tool versions
RUN set -eu; \
    check_version() { \
        tool="$1"; \
        expected="$2"; \
        actual="$3"; \
        if [ "$actual" != "$expected" ]; then \
            printf 'Version mismatch for %s: expected %s, installed %s\n' \
                "$tool" "$expected" "$actual" >&2; \
            exit 1; \
        fi; \
        printf 'Verified %s=%s\n' "$tool" "$actual"; \
    }; \
    check_version bcftools "${BCFTOOLS_VERSION}" \
        "$(bcftools --version | awk 'NR==1 {print $2}')"; \
    check_version tabix "${TABIX_VERSION}" \
        "$(tabix --version | awk 'NR==1 {print $3}')"; \
    check_version bgzip "${TABIX_VERSION}" \
        "$(bgzip --version | awk 'NR==1 {print $3}')"; \
    check_version vcfanno "${VCFANNO_VERSION}" \
        "$(vcfanno -h 2>&1 | awk '/^vcfanno version/ {print $3}')"

# Validate installed tools and Python dependencies
RUN python --version && \
    python -m pip check && \
    python -c "import numpy, Bio, pandas, pysam, sklearn, scipy, vcf2pandas, matplotlib"

LABEL org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.source="https://github.com/CCICB/introme" \
      org.ccicb.introme.runtime="cpu" \
      org.ccicb.introme.python.version="${PYTHON_VERSION}" \
      org.ccicb.introme.bcftools.version="${BCFTOOLS_VERSION}" \
      org.ccicb.introme.tabix.version="${TABIX_VERSION}" \
      org.ccicb.introme.vcfanno.version="${VCFANNO_VERSION}" \
      org.ccicb.introme.numpy.version="${NUMPY_VERSION}" \
      org.ccicb.introme.biopython.version="${BIOPYTHON_VERSION}" \
      org.ccicb.introme.pandas.version="${PANDAS_VERSION}" \
      org.ccicb.introme.pysam.version="${PYSAM_VERSION}" \
      org.ccicb.introme.scikit-learn.version="${SCIKIT_LEARN_VERSION}" \
      org.ccicb.introme.scipy.version="${SCIPY_VERSION}" \
      org.ccicb.introme.vcf2pandas.version="${VCF2PANDAS_VERSION}" \
      org.ccicb.introme.matplotlib.version="${MATPLOTLIB_VERSION}"
