## This repository contains code to process Single Cell RNA sequencing data generated using the QIAseqUltraplexRNA kit. 
```bash
# get the docker image with all the dependencies
sudo docker pull rpadmanabhan9/qiaseq-sc-rna:dev

# create a container with a shell mounting appropriate directories. please mount the directory with your fastq files here. the output files will be created in the same directory
sudo docker run -i -v /path/to/your_fastq_loc/:/path/to/your_fastq_loc/ -v /path/to/genome_dir/:/path/to/genome_dir/ rpadmanabhan9/qiaseq-sc-rna:dev bash

# clone this repository
cd /srv/qgen/code/
git https://github.com/qiaseq/qiaseq-singlecell-rna.git

# Run this pipeline as follows : 
cd qiaseq-singlecell-rna
EXPORT_LUIGI_CONFIG_PATH='pipeline.cfg'
luigid --background &
# for single cell experiment
PYTHONPATH="" luigi --module single_cell_rnaseq ClusteringAnalysis --samples-cfg samples.cfg --workers 22
# for low-input experiment
PYTHONPATH="" luigi --module single_cell_rnaseq WriteExcelSheet --samples-cfg samples.cfg --workers 22
```

### A few important parameters which you have to update in the pipeline.cfg file :
```
genome_dir : path to STAR genome index which is mounted on the container, i.e. /path/to/genome_dir/
primer_file : path to the primer file corresponding to your catalog number for a targeted sequencing experiment
annotation_gtf : path gencode annotation file , same as the one used to build the genome index with STAR
ercc_bed : path to the ERCC spike-in primer file
cell_index_file : path to the cell indices used in your experiment. File should have the oligo sequence of the indices, 1 per line

Please update any other parameters in the [DEFAULT] and [config] section as it applies to your sequencing experiment.

You can find the cell_index_file , primer_file , ercc_bed in the rpadmanabhan9/qiaseq-sc-rna:dev Docker Image under this path :
/home/qiauser/pipeline_data/

We have hosted some STAR indices in this gcp bucket : qiaseq-star-indices. To mount this bucket to your local system using gcsfuse :
gcsfuse qiaseq-star-indices /path/to/genome_dir/

You can find the annotation_gtf inside each organism's respective directory in the bucket. 
```
Update ***samples.cfg*** as it applies to your fastq files. 

Description of the analysis pipeline and output files can be found in the ***QIAseqUltraplexRNA_README.pdf file***.

We recommend you run this on a machine with atleast  60 GB RAM and 16 cores. To increase/decrease parallelism please tweak the **workers** command line option and/or the **num_cores** parameter in pipeline.cfg

