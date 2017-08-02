## This repository contains code to process Single Cell RNA sequencing data

### Run this pipeline as follows : 
PYTHONPATH="" EXPORT_LUIGI_CONFIG_PATH='./luigi.cfg' luigi --module single_cell_rnaseq JoinCountFiles --workers 8
