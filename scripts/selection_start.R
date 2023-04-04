#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(readr))
quiet(library(dplyr))
quiet(library(tidyr))
quiet(library(tibble))
quiet(library(stringr))
quiet(library(ggplot2))
quiet(library(RcppCNPy))

# the min frequency threshold to find
DELTA_MAF <- 0.1

clues_trajectory <- function(model_path) {

  # load the model data
  epochs <- npyLoad(paste0(model_path, ".epochs.npy"))
  freqs <- npyLoad(paste0(model_path, ".freqs.npy"))
  logpost <- npyLoad(paste0(model_path, ".post.npy"), dotranspose = F)

  # add column names
  colnames(logpost) <- paste0("V", seq(ncol(logpost)))

  model <- as_tibble(logpost) %>%

    # convert posterior densities from log-likelihoods
    exp() %>%

    # add the frequency labels
    add_column(freq = freqs, .before = 2) %>%

    # add the title heights (and a little padding to make sure there are no gaps)
    # tile heights are not equal because there is higher sampling density near 0 and 1
    add_column(height = diff(c(0, freqs)) + 1e-4, .before = 3) %>%

    # pivot the columns into long format
    pivot_longer(-c(freq, height), names_to = "epoch", values_to = "density") %>%

    # convert the column names into epochs (and switch the direction of time)
    mutate(epoch = -epochs[as.numeric(gsub("V", "", epoch))]) %>%

    # sort by epoch age
    arrange(epoch)

  model
}

# load the LCT sweep report
sweep <- read_tsv("lct-sweep/lct-sweep-snps-report.tsv", guess_max = 1500, col_types = cols()) %>%
    filter(p.value <= 5e-8) %>%
    rename(rsid=rsid_x)

# find the earliest epoch in which the frequency is above 10% with >0.5 density
models <- list()
for (i in 1:nrow(sweep)) {
    prefix <- paste0("clues/ancestral_paths_v3/all/", sweep[i, ]$rsid, "/ancestral_paths_v3-all-", sweep[i, ]$rsid, "-ancient-ALL-any")

    # load the CLUES model, and extract the maximum likelihood path
    model <- clues_trajectory(prefix) %>%
        group_by(epoch) %>%
        top_n(1, density) %>%
        ungroup() %>%
        arrange(epoch)

    starting_freq <- model %>% group_by() %>% slice_min(epoch) %>% pull(freq)

    if (sweep[i, ]$s > 0) {
        model <- model %>% filter(freq >= (starting_freq + DELTA_MAF))
    } else {
        # handle negatively selected alleles
        model <- model %>% filter(freq <= (starting_freq - DELTA_MAF))
    }

    model <- model %>%
        slice(which.min(epoch)) %>%
        mutate(rsid = sweep[i, ]$rsid, age_maf10 = -epoch * 28) %>%
        select(rsid, epoch_maf10=epoch, age_maf10)

    models <- append(models, list(model))
}

# bind_rows(models) %>% arrange(desc(age_maf10)) %>% mutate(rank=row_number()) %>% filter(rsid %in% c("rs4988235", "rs1438307")) %>% View()

sweep <- sweep %>%
    inner_join(bind_rows(models)) %>%
    arrange(desc(age_maf10))

write_tsv(sweep, "lct-sweep-snps-report-ages.tsv")
