#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

# code adapted from https://www.r-graph-gallery.com/101_Manhattan_plot.html

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(bedr))
quiet(library(data.table))
quiet(library(dplyr))
quiet(library(ggplot2))
quiet(library(ggrastr))
quiet(library(ggrepel))
quiet(library(stringr))

# load some helper functions
source("scripts/clues_utils.R")

# get the command line arguments
p <- arg_parser("Print a Manhattan plot of the CLUES data")
p <- add_argument(p, "--data", help = "CLUES data")
p <- add_argument(p, "--pairs", help = "GWAS and Control SNP pairings")
p <- add_argument(p, "--known", help = "List of known SNPs under selection")
p <- add_argument(p, "--unmapped", help = "List of SNPs unmapped by Relate")
p <- add_argument(p, "--flipped", help = "List of SNPs flipped by Relate")
p <- add_argument(p, "--output", help = "Filename of the PDF to render")
p <- add_argument(p, "--facet-x", help = "The x-axis facet")
p <- add_argument(p, "--facet-y", help = "The y-axis facet")

argv <- parse_args(p)

# use a genome-wide significance threshold
p.signif <- 5e-8
p.suggest <- 1e-5

# load all the data
data <- fread(argv$data, header = T, sep = "\t")
pairs <- fread(argv$pairs, header = T, sep = "\t")
unmapped <- fread(argv$unmapped, header = T, sep = "\t")
flipped <- fread(argv$flipped, header = T, sep = "\t")

# number of labeled significant SNPs per chrom
num_label <- 1

# determine if SNPs are GWAS or controls
data$type <- ifelse(data$rsid %in% pairs$gwas, "gwas", "control")

# set facet ordering
data$mode <- factor(data$mode, levels = c("modern", "ancient"), labels = c("Modern", "Ancient"))
data$type <- factor(data$type, levels = c("gwas", "control"), labels = c("GWAS", "Control"))
data$ancestry <- factor(data$ancestry, levels = c("ALL", "WHG", "EHG", "CHG", "ANA"))

# drop any modern SNPs that were flipped or unmapped
data <- data %>%
  mutate(pos = paste0(chrom, ":", start)) %>%
  filter(!(mode == "Modern" & (pos %in% unmapped$pos | pos %in% flipped$pos)))

# filter the data to only include SNPs that are in the pairs list
data <- data[data$rsid %in% c(pairs$gwas, pairs$neutral), ]

# annotate SNPs that were found in another study
data <- annotate_known_snps(data, argv$known, p.signif)

facets <- c(argv$facet_y, argv$facet_x)

# filter the dataset based on the facet columns
if (!("mode" %in% facets)) {
  data <- data %>%
    filter(mode == "Ancient")
}

if (!("ancestry" %in% facets)) {
  data <- data %>%
    filter(ancestry == "ALL")
}

if (!("type" %in% facets)) {
  data <- data %>%
    filter(type == "GWAS")
}

if ("Modern" %in% data$mode) {
  # drop any SNPs where we don't have a completed `modern` GWAS model (because Relate does not map all mutations)
  moderns <- data[data$mode == "Modern", ]$rsid
  pairs <- pairs[pairs$gwas %in% moderns & pairs$neutral %in% moderns, ]
}

# make the plot
plt <- manhattan_plot(data, facets, p.signif, p.suggest)

# png(file = argv$output, width = 16*1.3, height = 9*1.3, units = 'in', res = 300)
pdf(file = argv$output, width = 16 * 1.3, height = 9 * 1.3)
plt
dev <- dev.off()
