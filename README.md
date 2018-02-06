## This repository contains code to process Single Cell RNA sequencing data generated using the QIAseqUltraplexRNA kit. 
```bash
# Run this pipeline as follows : 
EXPORT_LUIGI_CONFIG_PATH='pipeline.cfg'
PYTHONPATH="" luigi --module single_cell_rnaseq ClusteringAnalysis --samples-cfg samples.cfg --workers 22
```


Description of the analysis pipeline and output files can be found in the QIAseqUltraplexRNA_README.pdf file.

You can use this Docker image to get an environment with all the dependencies :  ** rpadmanabhan9/qiaseq-sc-rna **

We recommend you run this on a machine with atleast  60 GB RAM and 16 cores. Please tweak the workers parameter appropriately. 
