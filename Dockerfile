
FROM biocontainers/biocontainers:latest

LABEL maintainer Raghavendra Padmanabahn <raghavendra.padmanabhan@qiagen.com>

## Create directory structure ##
USER root
RUN mkdir -p /srv/qgen/code/ && \
    mkdir -p /srv/qgen/data/ && \
    mkdir -p /srv/qgen/example/
    mkdir -p ~/luigi_scheduler_new/

## Install programs and libraries ##
RUN conda install star
RUN conda install pysam
RUN pip install pathos
RUN pip install intervaltree
RUN pip install regex
RUN pip install editdistance
RUN pip install biopython
RUN pip install guppy
RUN pip install luigi

## Add data ##
ADD scrna-seq/Primers_bed/ /srv/qgen/data/Primers_bed/
ADD scrna-seq/Cell_Index/ /srv/qgen/data/Cell_Index/
ADD scrna-seq/gencode.v23.basic.annotation.noPound.gtf.gz /srv/qgen/data/
ADD scrna-seq/human38_gencode/ /srv/qgen/data/human38_gencode/

## Start luigi ##
CMD luigid --background --pidfile ~/luigi_scheduler_new/luigi_pid --logdir ~/luigi_scheduler_new/ --state-path ~/luigi_scheduler_new/luigi_state
#EXPOSE 8082

## Environment Variables ##
ENV PYTHONPATH $PYTHONPATH:/opt/conda/lib/python2.7/site-packages/:
