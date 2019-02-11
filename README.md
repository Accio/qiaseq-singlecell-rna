Process Single Cell RNA sequencing data generated using the QIAseqUltraplexRNA kit
===
Jitao David Zhang, Feb 2019
:
Adapted from [https://github.com/qiaseq/qiaseq-singlecell-rna]. The adapted version works with Python3 + conda, without dependency on Docker.


## Major edits in the fork

* The code now works within an Conda environment with local intallation of STAR. No dependency on Docker is needed.
* intervaltree version 3.x is now used, use '.overlap' instead of '.search'
* Mixed indentation errors are fixed from the python code
* No shared memory is used given [apparent difficulty and potential issues](https://github.com/alexdobin/STAR/issues/277)

## Run the pipeline

```bash
sbatch qiagen-forSlurm.bash
```

## Appendix

## intervaltrees .search: .envelop or .overlap

Using .envelop caused many fewer reads annotated than using .overlap in `find_genes.py`. In fact, when I checked the `intervaltree` package in the docker file at [https://github.com/qiaseq/qiaseq-singlecell-rna](https://github.com/qiaseq/qiaseq-singlecell-rna), it turned out the `search` function had `strict=False`, suggesting that `overlap` should be used rather than `envelop`.
