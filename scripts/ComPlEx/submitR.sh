#!/bin/bash
set -eux

## variables
proj=u2020001
mail=elena.vanzalen@umu.se
rscript=scripts/ComPlEx/ComPlEx.R
out=data/ComPlEx

## create the out dir
if [ ! -d $out ]; then
    mkdir -p $out
fi

## execute
fnam=ComPlEx_C
sbatch -A $proj -t 2-00:00:00 --mail-user=$mail -e $out/$fnam.err -o $out/$fnam.out \
  -J $fnam -n 4 RunR.sh $rscript