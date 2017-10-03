
FROM biocontainers/biocontainers:latest

LABEL maintainer Raghavendra Padmanabahn <raghavendra.padmanabhan@qiagen.com>

## Create directory structure ##
USER root
RUN apt-get -y update && \
    apt-get -y install r-base

RUN mkdir -p /srv/qgen/code/
    
## Python dependencies
RUN conda install star
RUN pip install pysam
RUN pip install pathos
RUN pip install intervaltree
RUN pip install regex
RUN pip install editdistance
RUN pip install biopython
RUN pip install guppy
RUN pip install luigi
RUN pip install pandas

## R dependencies
#setup R configs
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('BASiCS')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('Rtsne')"
RUN Rscript -e "install.packages('MAST')"
RUN Rscript -e "install.packages('scde')"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('ggrepel')"
RUN Rscript -e "install.packages('cluster')"

## Environment Variables ##
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:
