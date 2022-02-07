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
quiet(library(tibble))
quiet(library(RcppCNPy))

# load some helper functions
source("scripts/clues_utils.R")

# load the clues report
report <- read_tsv("clues/ancestral_paths_new-all-clues_report.tsv", col_types = cols()) %>%
  rename(rsid = rsid_x) %>%
  select(c("rsid", "mode", "ancestry")) %>%
  drop_na()

freq <- bind_rows(
  apply(report, 1, function(row) {
    # extract the maximum posterior final frequency
    clues_load_data("ancestral_paths_new", "all", row["rsid"], row["mode"], row["ancestry"]) %>%
      filter(epoch == 0) %>%
      group_by(rsid, epoch) %>%
      top_n(1, density) %>%
      ungroup() %>%
      rename(freq_final = freq) %>%
      mutate(mode = row["mode"], ancestry = row["ancestry"]) %>%
      select(rsid, mode, ancestry, freq_final)
  })
)

write_tsv(freq, "clues/ancestral_paths_new-all-clues_report-freq_final.tsv")
