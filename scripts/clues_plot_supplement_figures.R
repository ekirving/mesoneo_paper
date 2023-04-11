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
p <- arg_parser("Plot the CLUES supplement figures")
p <- add_argument(p, "--data", help = "CLUES report", default = "clues/ancestral_paths_v3-all-filtered-clues_report.tsv")
p <- add_argument(p, "--pairs", help = "GWAS and neutral SNP pairings", default = "variants/ancestral_paths_v3-all-pairs.tsv")
p <- add_argument(p, "--unmapped", help = "List of modern SNPs unmappable by Relate", default = "relate/1000G_phase3-FIN_GBR_TSI-popsize-allsnps_unmapped.tsv.gz")
p <- add_argument(p, "--flipped", help = "List of modern SNPs flipped by Relate", default = "relate/1000G_phase3-FIN_GBR_TSI-popsize-allsnps_flipped.tsv.gz")
p <- add_argument(p, "--peaks-mod", help = "CLUES selection peaks from the modern analysis", default = "clues/ancestral_paths_v3-all-modern-ALL-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-all", help = "CLUES selection peaks from the pan-ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-ALL-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-whg", help = "CLUES selection peaks from the WHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-WHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-ehg", help = "CLUES selection peaks from the EHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-EHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-chg", help = "CLUES selection peaks from the CHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-CHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-ana", help = "CLUES selection peaks from the ANA ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-ANA-filtered-sweeps.tsv")
p <- add_argument(p, "--mathieson", help = "Selection peaks from Mathieson et al. 2015", default = "mathieson/mathieson-sweeps.tsv")
p <- add_argument(p, "--output", help = "Output prefix", default = "figs/supplement/ancestral_paths_v3-all")

argv <- parse_args(p)

# load some helper functions
source("scripts/clues_utils.R")

# load all the CLUES results
clues <- load_clues(argv$data, argv$pairs, argv$unmapped, argv$flipped)

# load all the CLUES peaks
modern_all_peaks <- read_tsv(argv$peaks_mod, col_types = cols())
ancient_all_peaks <- read_tsv(argv$peaks_all, col_types = cols())
ancient_whg_peaks <- read_tsv(argv$peaks_whg, col_types = cols())
ancient_ehg_peaks <- read_tsv(argv$peaks_ehg, col_types = cols())
ancient_chg_peaks <- read_tsv(argv$peaks_chg, col_types = cols())
ancient_ana_peaks <- read_tsv(argv$peaks_ana, col_types = cols())

# peaks for the reported p-values in Mathieson_et_al_2015
mathieson_peaks <- read_tsv(argv$mathieson, col_types = cols())

# color brewer https://colorbrewer2.org/#type=qualitative&scheme=Set2&n=5
ancestry_colors <- c(
  "ALL" = "#66c2a5",
  "WHG" = "#fc8d62",
  "EHG" = "#8da0cb",
  "CHG" = "#e78ac3",
  "ANA" = "#a6d854"
)

# get all the modern peaks
peaks_mod <- annotate_peaks_by_type(filter(clues, mode == "Modern"), modern_all_peaks, mathieson_peaks)

# get all the ancient ancestry models
ancestral <- filter(clues, mode == "Ancient")

# annotate ancestry specific peaks
peaks_all <- annotate_peaks_by_type(filter(ancestral, ancestry == "ALL"), ancient_all_peaks, mathieson_peaks)
peaks_whg <- annotate_peaks_by_type(filter(ancestral, ancestry == "WHG"), ancient_whg_peaks, mathieson_peaks)
peaks_ehg <- annotate_peaks_by_type(filter(ancestral, ancestry == "EHG"), ancient_ehg_peaks, mathieson_peaks)
peaks_chg <- annotate_peaks_by_type(filter(ancestral, ancestry == "CHG"), ancient_chg_peaks, mathieson_peaks)
peaks_ana <- annotate_peaks_by_type(filter(ancestral, ancestry == "ANA"), ancient_ana_peaks, mathieson_peaks)

# only show the marginal peaks in these plots
peaks_mod$data <- peaks_mod$data %>% mutate(merged_peaks = ifelse(peak != "absent", region, NA))
peaks_all$data <- peaks_all$data %>% mutate(merged_peaks = ifelse(peak != "absent", region, NA))
peaks_whg$data <- peaks_whg$data %>% mutate(merged_peaks = ifelse(peak != "absent", region, NA))
peaks_ehg$data <- peaks_ehg$data %>% mutate(merged_peaks = ifelse(peak != "absent", region, NA))
peaks_chg$data <- peaks_chg$data %>% mutate(merged_peaks = ifelse(peak != "absent", region, NA))
peaks_ana$data <- peaks_ana$data %>% mutate(merged_peaks = ifelse(peak != "absent", region, NA))

# merge all the ancestry stratified data back together
ancestral <- bind_rows(
    peaks_all$data,
    peaks_whg$data,
    peaks_ehg$data,
    peaks_chg$data,
    peaks_ana$data
)

# use a genome-wide significance threshold
p.genomewide <- 5e-8

# ------------------------------------------------------------------------------------------------------------------------------------------------------
# Modern 1000G data
# ------------------------------------------------------------------------------------------------------------------------------------------------------

plt_modern <- manhattan_plot(
    peaks_mod$data, 
    c("type"), 
    p.genomewide, 
    wrap = 1,
    colors = ancestry_colors, 
    gene_names = NA,
    size_snps = FALSE
) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("CLUES Modern analysis")

# plt_modern

ggsave(filename = paste0(argv$output, "-modern.png"), plot = plt_modern, width = 16, height = 8, units = "cm", scale = 1.5)

# ------------------------------------------------------------------------------------------------------------------------------------------------------
# Pan-ancestry aDNA data
# ------------------------------------------------------------------------------------------------------------------------------------------------------

plt_ancient <- manhattan_plot(
    peaks_all$data, 
    c("type"), 
    p.genomewide, 
    wrap = 1,
    colors = ancestry_colors, 
    gene_names = NA,
    size_snps = FALSE
) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("CLUES aDNA analysis")

# plt_ancient

ggsave(filename = paste0(argv$output, "-ancient.png"), plot = plt_ancient, width = 16, height = 8, units = "cm", scale = 1.5)

# ------------------------------------------------------------------------------------------------------------------------------------------------------
# Ancestry stratified aDNA data
# ------------------------------------------------------------------------------------------------------------------------------------------------------

max_height <- -log10(min(ancestral$p.value))

plt_gwas <- manhattan_plot(
    filter(ancestral, type=="GWAS"), 
    c("ancestry", "type"), 
    p.genomewide, 
    wrap = 1,
    colors = ancestry_colors, 
    gene_names = NA,
    size_snps = FALSE
) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("CLUES ancestral path GWAS SNPs") +
    expand_limits(y = max_height)

# plt_gwas

ggsave(filename = paste0(argv$output, "-ancestral-gwas.png"), plot = plt_gwas, width = 16, height = 20, units = "cm", scale = 1.5)

plt_control <- manhattan_plot(
    filter(ancestral, type=="Control"), 
    c("ancestry", "type"), 
    p.genomewide, 
    wrap = 1,
    colors = ancestry_colors, 
    gene_names = NA,
    size_snps = FALSE
) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ggtitle("CLUES ancestral path Control SNPs") +
    expand_limits(y = max_height)

# plt_control

ggsave(filename = paste0(argv$output, "-ancestral-control.png"), plot = plt_control, width = 16, height = 20, units = "cm", scale = 1.5)

# ------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------


