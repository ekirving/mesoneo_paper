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

# load the inversion calls from Andres
inv_calls <- read_delim("data/andres/inversion.calls.andres.txt", " ")

# load the sample metadata
meta <- read_tsv("data/ancestral_paths_v3/ancestral_paths_merged_filtered_age.sampleInfo.tsv")

# convert the inversion calls into CLUES input format
inv_calls %>%
  drop_na() %>%
  inner_join(meta, by = c("sample" = "sampleId")) %>%
  arrange(age) %>%
  mutate(call = case_when(
    h2_gt == 0 ~ "0.000000 -inf -inf",
    h2_gt == 1 ~ "-inf 0.000000 -inf",
    h2_gt == 2 ~ "-inf -inf 0.000000"
  )) %>%
  mutate(clues = sprintf("%.6f %s", age / 28, call)) %>%
  select(clues) %>%
  write_delim("clues/andres/inv17/inv17-clues.ancient", col_names = F, quote_escape = F, delim = "")
