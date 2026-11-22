FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

ARG IMAGE_VERSION=1.0.1

# Full APT package versions, including Ubuntu/Debian revisions.
# Defaults target Ubuntu 22.04 amd64.
ARG BCFTOOLS_VERSION=1.13-1
ARG BEDTOOLS_VERSION=2.30.0+dfsg-2ubuntu0.1
ARG TABIX_VERSION=1.13+ds-2build1
ARG VCFANNO_VERSION=0.3.3+ds-1ubuntu0.1

# Bioinformatics tools and runtime shell utilities
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        "bcftools=${BCFTOOLS_VERSION}" \
        "bedtools=${BEDTOOLS_VERSION}" \
        "tabix=${TABIX_VERSION}" \
        "vcfanno=${VCFANNO_VERSION}" \
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

# Verify installed package versions match the arguments used for labels.
# Single quotes preserve ${Version} as a dpkg-query format expression.
RUN set -eu; \
    for spec in \
        "bcftools=${BCFTOOLS_VERSION}" \
        "bedtools=${BEDTOOLS_VERSION}" \
        "tabix=${TABIX_VERSION}" \
        "vcfanno=${VCFANNO_VERSION}"; do \
        package="${spec%%=*}"; \
        expected="${spec#*=}"; \
        actual="$(dpkg-query -W -f='${Version}' "$package")"; \
        if [ "$actual" != "$expected" ]; then \
            printf 'Version mismatch for %s: expected %s, installed %s\n' \
                "$package" "$expected" "$actual" >&2; \
            exit 1; \
        fi; \
        printf 'Verified %s=%s\n' "$package" "$actual"; \
    done

CMD ["/bin/bash"]

LABEL org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.source="https://github.com/CCICB/introme" \
      org.ccicb.introme.runtime="cpu" \
      org.ccicb.introme.bcftools.version="${BCFTOOLS_VERSION}" \
      org.ccicb.introme.bedtools.version="${BEDTOOLS_VERSION}" \
      org.ccicb.introme.tabix.version="${TABIX_VERSION}" \
      org.ccicb.introme.vcfanno.version="${VCFANNO_VERSION}"

