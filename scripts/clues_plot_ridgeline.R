#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(dplyr))
quiet(library(purrr))
quiet(library(stringr))
quiet(library(jsonlite))
quiet(library(tidyr))
quiet(library(ggplot2))
quiet(library(ggridges))
quiet(library(viridis))
quiet(library(data.table))
quiet(library(scales))

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
report <- args[1]

# load all the models
data <- fread(report, header = T, sep = "\t")

# get a p-value for the test
data$p.value <- pchisq(2 * data$logLR, df = 1, lower.tail = FALSE)

# set a p-value cutoff
p_val <- 5e-8

# set s to zero for models above the p-value threshold
data[data$p.value > p_val]$s <- 0

# calculate absolute s
data$abs.s <- abs(data$s)

# set the type field
data$type <- ifelse(data$gwas == "", "Control", "GWAS")

# fill blank ancestry values
data$ancestry[data$ancestry == ""] <- "All"

# ancestry.order <- unique(data$ancestry)
# data$ancestry <- factor(data$ancestry, levels = ancestry.order)

###########################################################################################
set.seed(122)

# how many complete rsIDs are there
gwas_complete <- sample(filter(data, ancestry == "All" & mode == "modern" & type == "GWAS")$rsid)
gwas_num <- length(gwas_complete) # 1156
gwas_num

# trim the "control" SNPs to the same length
control_complete <- sample(filter(data, ancestry == "All" & mode == "modern" & type == "Control")$rsid)
control_num <- length(control_complete) # 1011
control_num

num_snps <- min(control_num, gwas_num)

if (control_num > gwas_num) {
  control_complete <- control_complete[1:gwas_num]
} else {
  gwas_complete <- gwas_complete[1:control_num]
}

# only get the complete models (inc. modern for all rsIDs)
data <- filter(data, rsid %in% c(gwas_complete, control_complete))

###########################################################################################
data2 <- data

data2$ancestry[data2$ancestry == "All"] <- str_to_sentence(data2$mode[data2$ancestry == "All"])

data2 <- data2 %>%
  pivot_wider(
    id_cols = c("rsid", "chrom", "pos", "genes", "anc/der", "epoch", "gwas", "type"),
    names_from = ancestry,
    values_from = abs.s
  ) %>%
  na.omit() %>%
  filter(type == "GWAS")

data2

# both
wilcox.test(data2$Ancient, data2$Modern, paired = TRUE)

###########################################################################################


###########################################################################################

all_snps <- data %>%
  filter(ancestry == "All") %>%
  group_by(mode) %>%
  mutate(mode_count = n()) %>%
  ungroup() %>%
  mutate(mode = paste0(mode, "; n=", mode_count))


# Ancient vs Modern
ggplot(filter(data, ancestry == "All"), aes(x = s, fill = mode)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~type, ncol = 1) +
  labs(title = paste0("Ancient vs Modern"), fill = "Mode")

ggplot(filter(data, ancestry == "All"), aes(x = s, fill = mode)) +
  geom_histogram(binwidth = 0.001, position = "identity", alpha = 0.5) +
  facet_wrap(~type, labeller = labeller(type = str_to_title), ncol = 1) +
  labs(title = paste0("GWAS vs Control | All SNPs"), fill = "SNP type")


# GWAS vs Control

# histogram
ggplot(all_snps, aes(x = s, fill = type)) +
  geom_histogram(binwidth = 0.001, position = "identity", alpha = 0.5) +
  facet_wrap(~mode, labeller = labeller(mode = str_to_title), ncol = 1) +
  labs(title = paste0("GWAS vs Control | All SNPs"), fill = "SNP type")

# density
ggplot(all_snps, aes(x = s, fill = type)) +
  geom_density(alpha = 0.25) +
  facet_wrap(~mode, labeller = labeller(mode = str_to_title), ncol = 1) +
  labs(title = paste0("GWAS vs Control | All SNPs"), fill = "SNP type")



selected_snps <- data %>%
  filter(ancestry == "All" & p.value < 0.05) %>%
  group_by(mode) %>%
  mutate(mode_count = n()) %>%
  ungroup() %>%
  mutate(mode = paste0(mode, "; n=", mode_count))

ggplot(selected_snps, aes(x = s, fill = type)) +
  geom_histogram(binwidth = 0.001, position = "identity", alpha = 0.5) +
  facet_wrap(~mode, ncol = 1) +
  labs(title = paste0("GWAS vs Control", " | Selected only"), fill = "SNP type")



###########################################################################################

# GWAS vs Control (Ancient + Ancestry paths)
# ggplot(filter(data, mode == 'ancient'), aes(x=s, fill=type)) +
#   geom_density(alpha=0.25) +
#   facet_wrap(~ancestry, ncol = 1) +
#   labs(title = paste0('GWAS vs Control (Ancient + Ancestry paths)', label_suffix), fill="SNP type")

all_snps <- data %>%
  filter(mode == "ancient") %>%
  group_by(ancestry) %>%
  mutate(ancestry_count = n()) %>%
  ungroup() %>%
  mutate(ancestry = paste0(ancestry, "; n=", ancestry_count))

ggplot(all_snps, aes(x = s, fill = type)) +
  geom_histogram(binwidth = 0.001, position = "identity", alpha = 0.5) +
  facet_wrap(~ancestry, ncol = 1) +
  labs(title = paste0("GWAS vs Control | Ancient only w/ Ancestry paths | All SNPs "), fill = "SNP type")

good_snps <- data %>%
  filter(mode == "ancient" & p.value < 0.05) %>%
  group_by(ancestry) %>%
  mutate(ancestry_count = n()) %>%
  ungroup() %>%
  mutate(ancestry = paste0(ancestry, "; n=", ancestry_count))

ggplot(good_snps, aes(x = s, fill = type)) +
  geom_histogram(binwidth = 0.001, position = "identity", alpha = 0.5) +
  facet_wrap(~ancestry, ncol = 1) +
  labs(title = paste0("GWAS vs Control | Ancient only w/ Ancestry paths | Selected only"), fill = "SNP type")


# # GWAS (Ancient + Ancestry paths)
# ggplot(filter(data, mode == 'ancient' & type == 'GWAS'), aes(x=s, fill=type)) +
#   geom_density(alpha=0.25) +
#   facet_wrap(~ancestry, ncol = 1) +
#   labs(title = paste0('GWAS SNPs (Ancient + Ancestry paths)', label_suffix), fill="SNP type")



###########################################################################################

# ggplot(data, aes(x = abs.s, y = ancestry)) +
#   geom_density_ridges() +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_discrete(expand = c(0, 0)) +
#   coord_cartesian(clip = "off") +
#   theme_minimal(base_size = 14) +
#   theme(axis.text.y = element_text(vjust = 0))


# # standard plot with all chains together
# ggplot(data, aes(x = abs.s, y = 'ancestry', fill = 0.5 - abs(0.5 - ..ecdf..))) +
#   stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
#                       jittered_points = TRUE, scale = 0.9, alpha = 0.7, color = 'grey',
#                       point_size = 0.4, point_alpha = 0.5, point_color = 'black',
#                       position = 'raincloud') +
#   # bandwidth = bandwidth
#   scale_fill_viridis(name = "Posterior", direction = -1) +
#   # theme_ridges(grid = FALSE) +
#   facet_wrap(~ancestry, ncol = 1) #+
#   # labs(title = toupper(measure)) +
#   # xlab(toupper(measure)) +
#   # limits
