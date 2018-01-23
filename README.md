## This repository contains code to process Single Cell RNA sequencing data
```bash
# Run this pipeline as follows : 
PYTHONPATH="" EXPORT_LUIGI_CONFIG_PATH='pipeline.cfg' luigi --module single_cell_rnaseq ClusteringAnalysis --samples-cfg samples.cfg --workers 22
