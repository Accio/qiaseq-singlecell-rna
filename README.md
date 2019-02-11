## This repository contains code to process Single Cell RNA sequencing data generated using the QIAseqUltraplexRNA kit. 


### Major edits in the fork

* The code now works within an Conda environment with local intallation of STAR. No dependency on Docker is needed.
* intervaltree version 3.x is now used
* Mixed indentation errors are fixed from the python code
* No shared memory is used given [apparent difficulty and potential issues](https://github.com/alexdobin/STAR/issues/277)

### QuickStart Example

The below commands will help you understand how to run our workflow. Please try to get this to run before attempting to modify any parameters or code. The example is from a low-input whole transcriptome experiment, the fastq files and config files are inside the docker image here : /home/qiauser/example/
```bash
# get the docker image with all the dependencies
sudo docker pull rpadmanabhan9/qiaseq-sc-rna:dev

# get genome index and annotation
# please note : we have stored the Genome Indices we use in the Google Cloud Bucket : qiaseq-star-indices
# link : https://console.cloud.google.com/storage/browser/qiaseq-star-indices
# you would need to install gsutil to download the files. Please see here for further details : https://cloud.google.com/storage/docs/gsutil_install#deb

gsutil -m cp -r gs://qiaseq-star-indices/human/STAR_index /path/to/your/genome_dir/

# create a container with a shell mounting appropriate directories.
# please specify /path/to/your_run_dir/ for storing the output files. This folder must be created by you.
# please specify /path/to/your/genome_dir/ . This is the same as the directory you downloaded the Genome index to

sudo docker run -it -v /path/to/your_run_dir/:/home/qiauser/rundir/ -v /path/to/your/genome_dir/:/home/qiauser/pipeline_data/STAR_index/ rpadmanabhan9/qiaseq-sc-rna:dev

# note : you are inside the container now
# clone this repository
cd /srv/qgen/code/
git clone https://github.com/qiaseq/qiaseq-singlecell-rna.git

# Run the pipeline, please modify --workers according the number of CPUs on your system. 
cd qiaseq-singlecell-rna
export LUIGI_CONFIG_PATH='/home/qiauser/example/pipeline.cfg'
luigid --background &
PYTHONPATH="" luigi --module single_cell_rnaseq WriteExcelSheet --samples-cfg /home/qiauser/example/samples.cfg --workers 22
```

### A few important parameters which you may have to update in the pipeline.cfg file according to your experiment :
```
seqtype           :    change to targetted for a targeted sequencing experiment
primer_file       :    path to the primer file corresponding to your catalog number for a targeted sequencing experiment
is_low_input      :    0 or 1 - whether this was a low-input experiment. 
                       Please change the WriteExcelSheet task in the luigi command line to ClusteringAnalysis 
                       for a singlecell experiment
annotation_gtf    :    path to gencode annotation file , same as the one used to build the genome index with STAR. Please check our
                       qiaseq-star-indices bucket for path to the annotation gtf file. 
cell_index_file   :    path to the cell indices used in your experiment. File should have the oligo sequence of the indices, 1 per line
                       path in docker image : /home/qiauser/pipeline_data/cell-index.list4.[96,384].txt
cell_indices_used :    comma delimeted list of cell indices used in your experiment - C1,C2,C3,C4,C5,.. 
                       to use all 96/384 to demux specify as : all
editdist          :    0 or 1 - editdistance mismatch for matching to known cell index set

Please update any other parameters in the [DEFAULT] and [config] section as it applies to your sequencing experiment.

```
Update ***samples.cfg*** as it applies to your fastq files. 

Description of the analysis pipeline and output files can be found in the ***QIAseqUltraplexRNA_README.pdf file***.

We recommend you run this on a machine with atleast  60 GB RAM and 16 cores. To increase/decrease parallelism please tweak the **workers** command line option and/or the **num_cores** parameter in pipeline.cfg

We are using luigi (https://github.com/spotify/luigi) to trigger our workflow. This allows for re-run on failure and fault tolerance among other features. The **single_cell_rnaseq.py** contains the task definitions for the various steps involved. The logic for the tasks themselves are under core/ in different modules.  

## intervaltrees .search: .envelop or .overlap

Using .envelop caused many fewer reads annotated than using .overlap in `find_genes.py`. In fact, when I checked the `intervaltree` package in the docker file at [https://github.com/qiaseq/qiaseq-singlecell-rna](https://github.com/qiaseq/qiaseq-singlecell-rna), it turned out the `search` function had `strict=False`, suggesting that `overlap` should be used rather than `envelop`.
