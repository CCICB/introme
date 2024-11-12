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
    bzip2

# Download and install HTSlib (which contains bgzip)
RUN curl -L https://github.com/samtools/htslib/releases/download/1.17/htslib-1.17.tar.bz2 -o htslib.tar.bz2 && \
    tar xvjf htslib.tar.bz2 && \
    cd htslib-1.17 && \
    ./configure && \
    make && \
    make install

# Clean up unnecessary files after installation
RUN rm -rf /var/lib/apt/lists/* htslib-1.17 htslib.tar.bz2

RUN pip install --upgrade pip
RUN pip install pysam Biopython pandas
