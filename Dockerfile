
FROM biocontainers/biocontainers:latest

LABEL maintainer Raghavendra Padmanabahn <raghavendra.padmanabhan@qiagen.com>

## Create directory structure ##
USER root
RUN apt-get -y update && \
    apt-get -y install r-base
    ## libraries for R packages
    apt-get install libnlopt-dev
    apt-get install libcairo2-dev
    apt-get install libgtk2.0-dev
    apt-get install xvfb
    apt-get install xauth
    apt-get install xfonts-base
    apt-get install libxt-dev
    apt-get install libssl-dev
    apt-get install libcurl4-openssl-dev    
    
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
RUN pip install natsort
RUN pip install xlsxwriter

## R dependencies
#setup R configs
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ggplot2')"
RUN Rscript -e "install.packages('dplyr')"
RUN Rscript -e "install.packages('Rtsne')"
RUN Rscript -e "install.packages('pheatmap')"

## Wierd version specific scde library installation
install_version("flexmix", version = "2.3-13", repos = "http://cran.us.r-project.org")
wget https://github.com/hms-dbmi/scde/archive/1.99.2.tar.gz
RUN Rscript -e "R CMD INSTALL 1.99.2.tar.gz"

## bioconductor related packages
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); biocLite('MAST');"
RUN Rscript -e "install.packages('stringr')"
RUN Rscript -e "install.packages('ggrepel')"
RUN Rscript -e "install.packages('cluster')"

## Some kwirks with this,currently copy from fdkbio07, have to install additional dependencies manually
RUN Rscript -e "install.packages('BASiCS')"

## Environment Variables ##
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:
