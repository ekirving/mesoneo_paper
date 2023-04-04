#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(dplyr))
quiet(library(ggplot2))
quiet(library(readr))
quiet(library(tidyr))

report <- read_tsv("variants/posterior_report_ancestral_paths_v3-neo_likelihoods.tsv") %>%

  # drop sites with no coverage one of the two datasets
  filter(epochs > 1) %>%

  # calculate the mean squared distances
  mutate(emd_sqr_mean = emd_tss / epochs, js_sqr_mean = js_tss / epochs)


report %>%
  select("rsid", "type", ends_with("_mean")) %>%
  pivot_longer(cols = -c("rsid", "type")) %>%

  # drop extreme values
  filter(value < 100) %>%

  # cap extreme values
  # mutate(value=pmin(value, 100)) %>%

  ggplot(aes(x = value, group = type, fill = type)) +
  geom_density(alpha = .4) +
  facet_wrap(~name, scales = "free", ncol = 2)


report %>%
  select("rsid", "type", ends_with("_mean")) %>%
  pivot_longer(cols = -c("rsid", "type")) %>%
  filter(name == "emd_mean") %>%
  ggplot(aes(x = value, group = type, fill = type)) +
  geom_density(alpha = .4) +
  geom_vline(xintercept = quantile(report$emd_mean, c(.95)), color = "orange") +
  geom_vline(xintercept = quantile(report$emd_mean, c(.99)), color = "red") +
  scale_x_log10()



quantile(report$emd_mean, c(.95)) # 10.06215
quantile(report$emd_mean, c(.99)) # 15.81759

quantile(report$emd_sqr_mean, c(.66)) # 24.22541
quantile(report$emd_sqr_mean, c(.95)) # 144.0587
quantile(report$emd_sqr_mean, c(.99)) # 364.0861


# ecdf(report$emd_mean)(25)
# ecdf(report$emd_sqr_mean)(25)

report %>%
  select("rsid", "type", ends_with("_mean")) %>%
  pivot_longer(cols = -c("rsid", "type")) %>%
  filter(name == "emd_sqr_mean") %>%
  filter(value > 144.0587) %>%
  View()
