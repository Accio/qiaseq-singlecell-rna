#!/bin/bash

#SBATCH -J qiagen
#SBATCH -o qiagen.out
#SBATCH -e qiagen.err
#SBATCH -c 4
#SBATCH --mem-per-cpu=40000
ml load Anaconda3/4.3.1
source activate QiagenSingleRna
ml load STAR
export LUIGI_CONFIG_PATH='/pstore/data/biomics/_pre_portfolio/_platform_evaluation/7788_LowInputEvaluation_3UPXQiagen/3UPXQiagenTest_19.12.2018/QiagenSinglecellRna/qiaseq-singlecell-rna/pipeline.cfg'
luigid &
PYTHONPATH='.' luigi --module single_cell_rnaseq WriteExcelSheet --samples-cfg /pstore/data/biomics/_pre_portfolio/_platform_evaluation/7788_LowInputEvaluation_3UPXQiagen/3UPXQiagenTest_19.12.2018/QiagenSinglecellRna/qiaseq-singlecell-rna/samples.cfg --workers 12
