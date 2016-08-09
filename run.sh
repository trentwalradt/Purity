#!/bin/bash -l
source ~/.bash_profile

module add R/3.1.3
module add samtools/1.1
module add gcc/4.9.2

Rscript $@
