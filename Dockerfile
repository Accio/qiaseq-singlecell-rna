
FROM biocontainers/biocontainers:latest

LABEL maintainer Raghavendra Padmanabahn <raghavendra.padmanabhan@qiagen.com>

## Create directory structure ##
USER root
RUN apt-get update
RUN mkdir -p /srv/qgen/code/
    
## Install programs and libraries ##
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

## Environment Variables ##
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:
