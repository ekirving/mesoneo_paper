#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# get the command line arguments
args=("$@")

datetime=$(date +%Y-%m-%d-%H%M)

# make a unique log file
logfile="nohup-${datetime}.out"

dataset=${args[0]:='ancestral_paths_v3'}
population=${args[1]:='all'}
start=${args[2]:=1}
stop=${args[3]:=22}
types=${args[4]:=gwas neutral}

# print the server name and start time to the log file
echo "SERVER: $HOSTNAME" >>${logfile}
echo "DATE: ${datetime}" >>${logfile}
echo "DATASET: ${dataset}-${population}" >>$logfile
echo "CHROMS: chr${start} - chr${stop}" >>$logfile

# load the conda environment
eval "$(conda shell.bash hook)"
conda activate mesoneo

MAX_ENSEMBL=15
MAX_MATPLOTLIB=20

if ! command -v free &> /dev/null; then
  # MacOS does not have the free command
  MAX_MEM=$(sysctl -a | awk '/^hw.memsize:/{print $2/(1024)^2}')
else
  # but linux does
  MAX_MEM=$(free -m | awk '/^Mem:/{print $2}')
fi

flags="--cores all "
flags+="--nolock "
flags+="--keep-going "
flags+="--printshellcmds "
flags+="--show-failed-logs "
flags+="--keep-incomplete "
flags+="--rerun-incomplete "
flags+="--reason "
#flags+="--restart-times 1 "
flags+="--resources mem_mb=${MAX_MEM} ensembl_api=${MAX_ENSEMBL} matplotlib=${MAX_MATPLOTLIB} "

for chr in $(seq ${start} ${stop}); do
  for type in $types; do
    for anc in ALL ANA CHG WHG EHG; do
      echo "Starting chr${chr} for ${anc} path(s) and ${type} SNPs..." >>$logfile
      (
        set -x
        snakemake ${flags} --config dataset=${dataset} population=${population} chr=${chr} ancestry=${anc} type=${type} -- ancestries
      ) &>>${logfile}
    done
  done
done

echo "DONE!" >>${logfile}
