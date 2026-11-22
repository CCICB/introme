FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Bioinformatics tools and runtime shell utilities
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        bcftools \
        bedtools \
        tabix \
        vcfanno \
        bash \
        gzip \
        coreutils \
        grep \
        sed \
        gawk \
        libc-bin \
    && rm -rf /var/lib/apt/lists/*

# Verify the expected executables exist
RUN set -eux; \
    for cmd in \
        bcftools bedtools tabix bgzip vcfanno \
        bash gunzip awk grep sort sed getconf; do \
        command -v "$cmd"; \
    done

CMD ["/bin/bash"]
