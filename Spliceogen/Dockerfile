FROM ubuntu:18.04

# Install OpenJDK-11
RUN apt-get update && \
    apt-get install -y openjdk-11-jre-headless && \
    apt-get clean;

RUN apt-get install bedtools

COPY . ./spliceogen
WORKDIR ./spliceogen

RUN chmod +x RUN.sh
ENTRYPOINT ["./RUN.sh"]
