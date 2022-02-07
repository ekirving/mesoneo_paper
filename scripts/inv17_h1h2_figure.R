#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(readr))
quiet(library(RcppCNPy))
quiet(library(tibble))
quiet(library(tidyr))
quiet(library(ggplot2))
quiet(library(directlabels))
quiet(library(ggpubr))
quiet(library(ggrepel))
quiet(library(shadowtext))
quiet(library(rjson))
quiet(library(dplyr))
quiet(library(zoo))

# load some helper functions
source("scripts/clues_utils.R")

dataset <- "ancestral_paths_new"
population <- "all"
mode <- "ancient"

ancestries <- c("ALL", "WHG", "EHG", "CHG", "ANA")

# load the all the chr17 inversion tag SNPs
tag_snps <- read_tsv("data/andres/inv17_h1h2_snps.txt")

# now fetch the inv17_h1h2 report, so we can plot the ones that were run
report <- read_tsv("clues/ancestral_paths_new-all-ancient-ALL-inv17_h1h2_report.tsv")

# join the dataframes
report <- report %>%
  left_join(tag_snps, by = c("rsid_x" = "SNP"))

# use the H1 allele to repolarize the trajectories
ancestral <- setNames(report$H1, as.character(report$rsid_x))

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# export the trajectory data

# load all the trajectories
traj <- bind_rows(
  lapply(ancestries, function(ancestry) {
    lapply(report$rsid_x, function(rsid) {
      clues_trajectory(dataset, population, rsid, mode, ancestry, smooth = 0, ancestral = ancestral[[rsid]]) %>% mutate(ancestry = ancestry)
    })
  })
)

# add the position back on
traj_report <- traj %>%
  inner_join(report %>% select(rsid_x, chrom, start), by = c("rsid" = "rsid_x")) %>%
  select(rsid, chrom, start, ancestry, epoch, freq, density)

write_tsv(traj_report, "inv17/ancestral_paths_new-all-ancient-inv17_h1h2-trajectories.tsv")

# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------


# make the per-ancesty CLUES plots
plt_clues <- lapply(ancestries, function(focal_ancestry) {
  clues_plot(dataset, population, report$rsid_x, mode, focal_ancestry, NA, geom = "line", title = focal_ancestry, ancestral = ancestral)
})

plt_inv <- ggarrange(plotlist = plt_clues, nrow = 1)

ggsave(filename = "inv17_h1h2_fig.png", plot = plt_inv, width = 16 * 10, heigh = 9 * 5, units = "mm", scale = 4, limitsize = FALSE)


# ------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------

# remove 1 problem SNP
report <- report %>%
  filter(!rsid_x %in% c("rs62058962"))

clues_imput <- clues_plot("ancestral_paths_new", population, report$rsid_x, mode, "ALL", NA, geom = "line", title = "Imputed", ancestral = ancestral)
# clues_imput

clues_likeli <- clues_plot("neo_likelihoods", population, report$rsid_x, mode, "ALL", NA, geom = "line", title = "Likelihood", ancestral = ancestral)
# clues_likeli

plt <- ggarrange(clues_imput, clues_likeli)
ggsave(filename = "inv17/inv17-imputed_vs_likelihood.png", plot = plt, width = 16 * 5 / 2, heigh = 9 * 5 / 2, units = "mm", scale = 4, limitsize = FALSE)

# load the EMD report
# emd <- read_tsv("variants/posterior_report_ancestral_paths_new-neo_likelihoods.tsv") %>%
#
#   # drop sites with no coverage one of the two datasets
#   filter(epochs > 1) %>%
#
#   # calculate the mean squared distances
#   mutate(emd_sqr_mean = emd_tss / epochs, js_sqr_mean = js_tss / epochs)
#
#
# report$rsid_x
