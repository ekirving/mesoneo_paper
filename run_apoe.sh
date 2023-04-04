#!/usr/bin/env bash

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT


# make the CLUES inputs
Rscript scripts/apoe_input.R

# the three different models to run
ISOFORMS=( APOE2 APOE3 APOE4 )

# run all the models
for isoform in "${ISOFORMS[@]}"; do

  python bin/clues/inference.py \
    --lik \
    --popFreq `cat "apoe/${isoform}.freq"` \
    --coal relate/1000G_phase3-FIN_GBR_TSI-popsize.coal \
    --ancientHaps "apoe/${isoform}.ancient" \
    --timeBins clues/ancestral_paths_v3-all-time.bins \
    --betaParam 0.5  \
    --out "apoe/${isoform}" &> "apoe/${isoform}.log"

  # extract the results
  python scripts/clues_parse_log.py \
    --rsid ${isoform} \
    --mode ancient \
    --ancestry ALL \
    --sex any \
    --log "apoe/${isoform}.log" \
    --out "apoe/${isoform}.json"

  num_samples=$(($(cat "apoe/${isoform}.ancient" | wc -l) / 2))

  # make the label
  echo "{\"title\":\"${isoform} (n=${num_samples})\", \"gwascat\": []}" \
    > "apoe/${isoform}-label.json"

  # plot the trajectory
  python scripts/clues_plot_trajectory.py \
    --gen-time 28 \
    --params "apoe/${isoform}.json" \
    --label  "apoe/${isoform}-label.json" \
    --ancestry ALL \
    --sex any \
    --ext png \
    "apoe/${isoform}" \
    "apoe/${isoform}" 2> /dev/null

done;

