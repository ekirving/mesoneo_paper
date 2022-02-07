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
quiet(library(ggrepel))
quiet(library(RcppCNPy))
quiet(library(tibble))
quiet(library(zoo))
quiet(library(directlabels))
quiet(library(ggpubr))
quiet(library(ggpubr))
quiet(library(rjson))

# load some helper functions
source("scripts/clues_utils.R")

p.genomewide <- 5e-8
p.bonferroni <- 1.35e-6

ancestries <- c("ALL", "WHG", "EHG", "CHG", "ANA")

# color brewer https://colorbrewer2.org/#type=qualitative&scheme=Set2&n=5
ancestry_colors <- c(
  "ALL" = "#66c2a5", # "grey"
  "WHG" = "#fc8d62",
  "EHG" = "#8da0cb",
  "CHG" = "#e78ac3",
  "ANA" = "#a6d854"
)

# color brewer https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=12
snp_colors <- c(
  "#a6cee3",
  "#1f78b4",
  "#b2df8a",
  "#33a02c",
  "#fb9a99",
  "#e31a1c",
  "#fdbf6f",
  "#ff7f00",
  "#cab2d6",
  "#6a3d9a",
  "#ffff99",
  "#b15928"
)

# fetch the inv17_h1h2 reports
region_data <- bind_rows(
  read_tsv("clues/ancestral_paths_new-all-ancient-ALL-inv17_h1h2_report.tsv"),
  read_tsv("clues/ancestral_paths_new-all-ancient-ANA-inv17_h1h2_report.tsv"),
  read_tsv("clues/ancestral_paths_new-all-ancient-CHG-inv17_h1h2_report.tsv"),
  read_tsv("clues/ancestral_paths_new-all-ancient-EHG-inv17_h1h2_report.tsv"),
  read_tsv("clues/ancestral_paths_new-all-ancient-WHG-inv17_h1h2_report.tsv"),
)

# fix `rsid` colname issue
region_data <- region_data %>%
  rename(rsid = rsid_x) %>%
  select(-rsid_y)

# add dummy entry for "merged_peaks"
region_data <- region_data %>%
  mutate(merged_peaks = "chr17:inv")

# get the top SNP from each of the marginal ancestries
top_marginal_snps <- region_data %>%
  group_by(ancestry) %>%
  slice_min(p.value, n = 3) %>%
  group_by(rsid) %>%
  summarise(n = length(rsid), p.sum = sum(-log10(p.value))) %>%
  arrange(desc(p.sum))

# get the list of top SNPs in the region
region_snps <- head(top_marginal_snps$rsid, n = length(snp_colors))

# pair the colours with the ordered top SNPs
label_colors <- setNames(snp_colors[1:length(region_snps)], region_snps)

# make region name human readable
region_name <- "chr17 inversion"

# determine the max height, to keep the y-axis in syc
max_height <- max(-log10(region_data$p.value))

# make the zoomed in manhattan plots
plt_col_zoom <- lapply(ancestries, function(focal_ancestry) {
  manhattan_plot(filter(region_data, ancestry == focal_ancestry), c("merged_peaks", "ancestry"), p.genomewide, p.bonferroni, top_snps = region_snps, show_strip = FALSE, wrap = 1, composite = FALSE, colors = ancestry_colors, label_colors = label_colors, region_name = region_name, max_height = max_height)
})

# load the all the chr17 inversion tag SNPs
tag_snps <- read_tsv("data/andres/inv17_h1h2_snps.txt")

# use the H1 allele to repolarize the trajectories
ancestral <- setNames(tag_snps$H1, as.character(tag_snps$SNP))

# make the per-ancesty CLUES plots
plt_col_clues <- lapply(ancestries, function(focal_ancestry) {

  # sort the SNPs by their marginal p-value
  ordered_snps <- region_data %>%
    filter(rsid %in% region_snps & ancestry == "ALL") %>%
    arrange(p.value) %>%
    pull(rsid)

  clues_plot("ancestral_paths_new", "all", ordered_snps, "ancient", focal_ancestry, label_colors, ancestral = ancestral)
})

plt_region <- ggarrange(
  ggarrange(plotlist = plt_col_zoom, nrow = 1, labels = ancestries, font.label = list(color = "grey45"), hjust = 0, vjust = 0),
  ggarrange(plotlist = plt_col_clues, nrow = 1),
  nrow = 2
) + theme(plot.title = element_text(hjust = 0.5)) + ggtitle(" ")

plt_region

ggsave(filename = "inv17/inv17_h1h2_fig-locus.png", plot = plt_region, width = 16 * 10, heigh = 9 * 5, units = "mm", scale = 4, limitsize = FALSE)
