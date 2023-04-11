#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2021, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(readr))
quiet(library(ggpubr))
quiet(library(gtools))

# # get the command line arguments
p <- arg_parser("Plot the CLUES simulated data figure")
p <- add_argument(p, "--data", help = "CLUES report", default = "clues/simulated_relate_painted-all-simulated-clues_report.tsv")
p <- add_argument(p, "--pairs", help = "GWAS and neutral SNP pairings", default = "variants/simulated_relate_painted-all-pairs.tsv")
p <- add_argument(p, "--peaks-all", help = "CLUES selection peaks from the pan-ancestry analysis", default = "clues/simulated_relate_painted-all-simulated-ALL-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-whg", help = "CLUES selection peaks from the WHG ancestry analysis", default = "clues/simulated_relate_painted-all-simulated-WHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-ehg", help = "CLUES selection peaks from the EHG ancestry analysis", default = "clues/simulated_relate_painted-all-simulated-EHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-chg", help = "CLUES selection peaks from the CHG ancestry analysis", default = "clues/simulated_relate_painted-all-simulated-CHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-ana", help = "CLUES selection peaks from the ANA ancestry analysis", default = "clues/simulated_relate_painted-all-simulated-ANA-filtered-sweeps.tsv")
p <- add_argument(p, "--merged", help = "Merged selection peaks across ancestries", default = "clues/simulated_relate_painted-all-simulated-filtered-sweeps-merged.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "figs/supplement/simulated_relate_painted-all-simulated.png")

argv <- parse_args(p)

# load some helper functions
source("scripts/clues_utils.R")

# load all the CLUES results
data <- read_tsv(argv$data, col_types = cols())
pairs <- read_tsv(argv$pairs, col_types = cols())

# load all the simulated peaks
simulated_all_peaks <- read_tsv(argv$peaks_all, col_types = cols())
simulated_whg_peaks <- read_tsv(argv$peaks_whg, col_types = cols())
simulated_ehg_peaks <- read_tsv(argv$peaks_ehg, col_types = cols())
simulated_chg_peaks <- read_tsv(argv$peaks_chg, col_types = cols())
simulated_ana_peaks <- read_tsv(argv$peaks_ana, col_types = cols())
simulated_merged_peaks <- read_tsv(argv$merged, col_types = cols())

# determine if SNPs are GWAS or controls
data$type <- ifelse(data$rsid %in% pairs$gwas, "gwas", "control")

# set facet ordering
data$mode <- factor(data$mode, levels = c("modern", "ancient"), labels = c("Modern", "Ancient"))
data$type <- factor(data$type, levels = c("gwas", "control"), labels = c("GWAS", "Control"))
data$ancestry <- factor(data$ancestry, levels = c("ALL", "WHG", "EHG", "CHG", "ANA"))

if (nrow(simulated_merged_peaks) == 0) {
  # there are no peaks to show
  data$region <- NA
  data$merged_peaks <- NA
}

# color brewer https://colorbrewer2.org/#type=qualitative&scheme=Set2&n=5
ancestry_colors <- c(
  "ALL" = "#66c2a5",
  "WHG" = "#fc8d62",
  "EHG" = "#8da0cb",
  "CHG" = "#e78ac3",
  "ANA" = "#a6d854"
)

# use a genome-wide significance threshold
p.genomewide <- 5e-8

# ------------------------------------------------------------------------------------------------------------------------------------------------------
# Ancestry stratified simulated data
# ------------------------------------------------------------------------------------------------------------------------------------------------------

plt_sim <- manhattan_plot(
    data, 
    c("ancestry", "type"), 
    p.genomewide, 
    wrap = 1,
    colors = ancestry_colors, 
    gene_names = NA,
    size_snps = FALSE
) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("CLUES ancestral path simulated SNPs")

# plt_sim

ggsave(filename = argv$output, plot = plt_sim, width = 16, height = 20, units = "cm", scale = 1.5)
