#!/bin/bash
#SBATCH -A hpc2n2026-283
#SBATCH -J complex_loo
#SBATCH -o logs/task_%a.out
#SBATCH -e logs/task_%a.out
ulimit -c 0
W=/proj/nobackup/hpc2nstor2025-104/complex_loo_20261007
cd $W
export PIXI_HOME=$W/.pixi_home PIXI_CACHE_DIR=$W/.pixi_cache
D=$(mktemp -d $W/runs/t${SLURM_ARRAY_TASK_ID}_XXXX)
cd $D && $W/pixi run --manifest-path $W/pixi.toml Rscript $W/task.R $W $SLURM_ARRAY_TASK_ID
rc=$?
cd $W && rm -rf $D
exit $rc
