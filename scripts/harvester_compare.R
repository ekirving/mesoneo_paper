#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(readr))
quiet(library(dplyr))
quiet(library(tidyr))
quiet(library(ggplot2))
quiet(library(bedr))

# load some helper functions
source("scripts/clues_utils.R")

ancestries <- c("ALL", "ANA", "CHG", "EHG", "WHG")

# load all the original peaks
peaks.orig <- havester_to_bed(bind_rows(lapply(ancestries, function(ancestry) {
  read_tsv(paste0("clues/ancestral_paths_v3-all-ancient-", ancestry, "-harvester.tsv"), col_types = cols()) %>%
    mutate(ancestry = ancestry)
})))

# load all the posterior_diff peaks
peaks.post <- havester_to_bed(bind_rows(lapply(ancestries, function(ancestry) {
  read_tsv(paste0("clues/ancestral_paths_v3-all-ancient-", ancestry, "-posterior_diff-harvester.tsv"), col_types = cols()) %>%
    mutate(ancestry = ancestry)
})))

# load all the refbias_Anc peaks
peaks.refbias_anc <- havester_to_bed(bind_rows(lapply(ancestries, function(ancestry) {
  read_tsv(paste0("clues/ancestral_paths_v3-all-ancient-", ancestry, "-refbias_Anc-harvester.tsv"), col_types = cols()) %>%
    mutate(ancestry = ancestry)
})))

# load all the refbias_Anc_1000g peaks
peaks.refbias_anc_1000g <- havester_to_bed(bind_rows(lapply(ancestries, function(ancestry) {
  read_tsv(paste0("clues/ancestral_paths_v3-all-ancient-", ancestry, "-refbias_Anc_1000g-harvester.tsv"), col_types = cols()) %>%
    mutate(ancestry = ancestry)
})))
