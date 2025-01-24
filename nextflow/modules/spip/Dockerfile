FROM --platform=linux/amd64 ubuntu

RUN apt-get update && \ 
	apt-get install -y wget && \
	DEBIAN_FRONTEND=noninteractive apt-get -qq install r-base r-base-dev

# Updates
RUN apt-get update
RUN apt-get install -y --no-install-recommends apt-utils
RUN apt-get -y install git

# Install the packages
RUN R -e "install.packages(c('foreach', 'doParallel', 'randomForest'), repos = 'http://cran.us.r-project.org')"

# Copy the repo
RUN git clone https://github.com/raphaelleman/SPiP
RUN cd ./SPiP/RefFiles \
  && wget -O transcriptome_hg38.RData https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg38.RData/download \
  && wget -O transcriptome_hg19.RData https://sourceforge.net/projects/splicing-prediction-pipeline/files/transcriptome/transcriptome_hg19.RData/download

WORKDIR /SPiP
