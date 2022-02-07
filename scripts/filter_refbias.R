#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(readr))
quiet(library(dplyr))
quiet(library(tidyr))

# get the command line arguments
p <- arg_parser("Filter the CLUES report based on Fj scores")
p <- add_argument(p, "--data", help = "CLUES report", default = "clues/ancestral_paths_new-all-clues_report.tsv")
p <- add_argument(p, "--gwas", help = "Fj scores for the GWAS SNPs", default = "refbias/Evan_GWAS_Anc_1000g.txt.gz")
p <- add_argument(p, "--neut", help = "Fj scores for the neutral SNPs", default = "refbias/Evan_NEUTRAL_Anc_1000g.txt.gz")
p <- add_argument(p, "--post", help = "Posterior frequency differences", default = "refbias/posterior-diff.tsv.gz")
p <- add_argument(p, "--output", help = "Output file")

argv <- parse_args(p)

# maximum absolute difference in inferred present-day frequency between models using 1000G populations and those using ancient DNA alone
POSTERIOR_MAX <- 0.1

# the range of Fj scores to exclude due to putative reference bias
Fj_MIN <- 0.5
Fj_MAX <- 1.0

# load all the data
data <- read_tsv(argv$data, col_types = cols())
gwas <- read_tsv(argv$gwas, col_types = cols())
neut <- read_tsv(argv$neut, col_types = cols())
post <- read_tsv(argv$post, col_types = cols())

cols <- colnames(data)

refbias <- bind_rows(gwas, neut)

# filter for Fj
data <- data %>%
  inner_join(refbias, by = c("chrom" = "Chr", "start" = "Pos")) %>%
  filter(Fj < Fj_MIN | Fj >= Fj_MAX) %>%
  select(all_of(cols))

# filter for posterior frequency
data <- post %>%
  filter(mode == "ancient" & ancestry == "ALL") %>%
  select(rsid, diff_abs, logLR_with_mod, logLR_no_mod) %>%
  inner_join(data, by = "rsid") %>%
  filter(diff_abs < POSTERIOR_MAX | logLR_no_mod > logLR_with_mod) %>%
  select(all_of(cols))

write_tsv(data, argv$output)
