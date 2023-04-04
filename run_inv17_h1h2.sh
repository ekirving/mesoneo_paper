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
batches=${args[2]:=100}

# print the server name and start time to the log file
echo "SERVER: $HOSTNAME" >>${logfile}
echo "DATE: ${datetime}" >>${logfile}
echo "DATASET: ${dataset}-${population}" >>$logfile

# load the conda environment
eval "$(conda shell.bash hook)"
conda activate mesoneo

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
flags+="--resources mem_mb=${MAX_MEM} "

for batch in $(seq 1 ${batches}); do
  echo "Starting batch ${batch} for inv17_h1h2..." >>$logfile
  (
    set -x
    snakemake ${flags} --batch inv17_h1h2=${batch}/${batches} --config dataset=${dataset} population=${population} -- inv17_h1h2
  ) &>>${logfile}
done

echo "DONE!" >>${logfile}
