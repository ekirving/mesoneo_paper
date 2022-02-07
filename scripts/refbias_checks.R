#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2021, University of Copenhagen
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

# load the CLUES reports
data.with_mod <- read_tsv("clues/ancestral_paths_new-all-clues_report.tsv", col_types = cols()) %>%
  filter(mode == "ancient") %>%
  select(c("rsid", "chrom", "start", "ancestry", "logLR", "s", "p.value", "info")) %>%
  mutate(log10p = -log10(p.value)) %>%
  mutate(vesion = "with_mod")

data.no_mod <- read_tsv("clues.no_mod/ancestral_paths_new-all-clues_report.tsv", col_types = cols()) %>%
  filter(mode == "ancient") %>%
  rename(rsid = rsid_x) %>%
  select(c("rsid", "chrom", "start", "ancestry", "logLR", "s", "p.value", "info")) %>%
  mutate(log10p = -log10(p.value)) %>%
  mutate(vesion = "no_mod")

# load all the refbias scores
gwas.anc <- read_tsv("refbias/Evan_GWAS_Anc.txt", col_types = cols())
gwas.mod <- read_tsv("refbias/Evan_GWAS_Anc_1000g.txt", col_types = cols())
neut.anc <- read_tsv("refbias/Evan_NEUTRAL_Anc.txt", col_types = cols())
neut.mod <- read_tsv("refbias/Evan_NEUTRAL_Anc_1000g.txt", col_types = cols())

# merge all the dataframes
refbias <- bind_rows(
  mutate(gwas.anc, type = "GWAS", subtype = "Anc"),
  mutate(gwas.mod, type = "GWAS", subtype = "Anc_1000g"),
  mutate(neut.anc, type = "NEUTRAL", subtype = "Anc"),
  mutate(neut.mod, type = "NEUTRAL", subtype = "Anc_1000g")
)

# what percentage of sites have Fj > 0.5
refbias %>%
  group_by(type, subtype) %>%
  summarize(count = sum(Fj > 0.5, na.rm = TRUE)) %>%
  mutate(percent = count / 36994)

# type      subtype    count  percent
# <chr>     <chr>      <int>  <dbl>
# 1 GWAS    Anc        1184   0.0320
# 2 GWAS    Anc_1000g  3320   0.0897
# 3 NEUTRAL Anc        1385   0.0374
# 4 NEUTRAL Anc_1000g  3459   0.0935

# the max values for Fj are very large, so thrshold it for making plots
thresh <- 2

refbias <- refbias %>%
  mutate(Fj_thres = ifelse(Fj > thresh, thresh, ifelse(Fj < -thresh, -thresh, Fj)))

# ------------------------------------------------------------------------------------------

# plot the distribution of the refbias scores
ggplot(data = refbias, aes(x = Fj_thres)) +
  facet_grid(subtype ~ type) +
  geom_density() +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "red")

ggsave("refbias/01-refbias-density.png", width = 16, height = 9)

# ------------------------------------------------------------------------------------------

# join the refbias scores to the p-values
data.refbias <- bind_rows(
  inner_join(data.with_mod, filter(refbias, subtype == "Anc_1000g"), by = c("chrom" = "Chr", "start" = "Pos")),
  inner_join(data.no_mod, filter(refbias, subtype == "Anc"), by = c("chrom" = "Chr", "start" = "Pos")),
)

# scatter plot of log10 p-value vs refbias scores
ggplot(data = data.refbias) +
  facet_grid(subtype ~ type) +
  geom_point(aes(x = Fj_thres, y = log10p)) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "red") +
  annotate("rect", xmin = 0.5, xmax = Inf, ymin = -log10(5e-8), ymax = Inf, alpha = 0.2, fill = "red")

ggsave("refbias/02-refbias-log10p-scatter.png", width = 16, height = 9)

# ------------------------------------------------------------------------------------------

# load the posterior frequencies of both sets of CLUES models (with modern genotype frequencies vs. without)
post.with_mod <- read_tsv("clues/ancestral_paths_new-all-clues_report-freq_final.tsv", col_types = cols())
post.no_mod <- read_tsv("clues.no_mod/ancestral_paths_new-all-clues_report-freq_final.tsv", col_types = cols())

data.both <- data.with_mod %>%
  inner_join(data.no_mod, by = c("rsid", "ancestry"), suffix = c("_with_mod", "_no_mod")) %>%
  select(rsid, ancestry, logLR_with_mod, logLR_no_mod)

# join the dfs
post <- post.with_mod %>%
  inner_join(post.no_mod, by = c("rsid", "mode", "ancestry"), suffix = c("_with_mod", "_no_mod")) %>%
  mutate(diff = freq_final_with_mod - freq_final_no_mod) %>%
  mutate(diff_abs = abs(diff)) %>%
  inner_join(data.both, by = c("rsid", "ancestry"))

# what percentage of SNPs are >10% off
post %>%
  group_by(mode, ancestry) %>%
  summarize(count = sum(diff_abs > 0.1, na.rm = TRUE)) %>%
  mutate(percent = count / 36994)

# mode      ancestry count percent
# <chr>     <chr>    <int>   <dbl>
# 1 ancient ALL        181 0.00489
# 2 ancient ANA        941 0.0254
# 3 ancient CHG       5131 0.139
# 4 ancient EHG       4223 0.114
# 5 ancient WHG       8646 0.234

write_tsv(post, "refbias/posterior-diff.tsv.gz")

# plot the distribution of frequency differences
ggplot(data = post) +
  facet_wrap(~ancestry) +
  geom_density(aes(x = diff), adjust = 3) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = "dashed", colour = "red") +
  scale_x_continuous(name = "Present-day frequency difference", breaks = seq(-1, 1, 0.2))

ggsave("refbias/03-posterior-diff.png", width = 16, height = 9)

# ------------------------------------------------------------------------------------------

# get a mapping table of rsIDs to chr:pos
rsids <- data.no_mod %>%
  select(rsid, chrom, start) %>%
  unique()

# join the refbias scores to the posterior-diff scores
post.refbias <- post %>%
  filter(ancestry == "ALL") %>%
  inner_join(rsids, by = "rsid") %>%
  inner_join(refbias, by = c("chrom" = "Chr", "start" = "Pos"))

# scatter plot of refbias vs posterior-diff
ggplot(data = post.refbias, aes(x = Fj_thres, y = diff)) +
  facet_wrap(~subtype) +
  geom_point() +
  geom_smooth(method = "lm", formula = "y ~ x") +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "red") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.1, ymax = Inf, alpha = 0.2, fill = "orange") +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = -0.1, alpha = 0.2, fill = "orange") +
  annotate("rect", xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2, fill = "red") +
  geom_label_repel(data = filter(post.refbias, Fj_thres > 0.5 & diff_abs > 0.1), mapping = aes(label = rsid), min.segment.length = 0, box.padding = 0.5, max.overlaps = 100) +
  geom_label_repel(data = filter(post.refbias, rsid == "rs80028338"), mapping = aes(label = rsid, color = "red"), min.segment.length = 0, box.padding = 0.5, max.overlaps = 100)

ggsave("refbias/04-posterior-diff_vs_refbias.png", width = 16, height = 9)

# ------------------------------------------------------------------------------------------
