# Use the latest Ubuntu LTS as the base image
FROM ubuntu:22.04

# Set environment variables to avoid interactive prompts during package installations
ENV DEBIAN_FRONTEND=noninteractive \
LANG=C.UTF-8 \
LC_ALL=C.UTF-8

# Install required system packages and clean up APT cache
RUN apt-get update && \
apt-get install -y --no-install-recommends openjdk-11-jre-headless bedtools git && \
rm -rf /var/lib/apt/lists/*

# Define build-time variables
ARG SPLICEOGEN_REPO=https://github.com/VCCRI/Spliceogen.git
ARG SPLICEOGEN_VERSION=a8a5d2d805781c7b25e4c788793e9fc19c94167a

# Clone the Spliceogen repository
RUN git init Spliceogen \
    && cd Spliceogen \
    && git remote add origin ${SPLICEOGEN_REPO} \
    && git fetch --depth 1 origin ${SPLICEOGEN_VERSION} \
    && git checkout FETCH_HEAD

# Set the working directory
WORKDIR /Spliceogen

# Make the RUN.sh script executable
RUN chmod +x RUN.sh

# Set the entrypoint to run the application
ENTRYPOINT ["./RUN.sh"]
