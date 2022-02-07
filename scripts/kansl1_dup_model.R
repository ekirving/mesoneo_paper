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

# load the duplication calls from Alma
dup_calls <- read_tsv("inv17/Ancient_17q_Inversion_Duplication_Pop_Age.txt", col_types = cols())

# load the inversion calls from Andres
inv_calls <- read_delim("data/andres/inversion.calls.andres.txt", " ", col_types = cols()) %>%

  # convert the calls into the format used by Alma
  mutate(andres_call = case_when(
    h2_gt == 0 ~ "H1H1",
    h2_gt == 1 ~ "H1H2",
    h2_gt == 2 ~ "H2H2"
  ))

# load the genotype calls for rs80028338 (the KANSL1 SNP under strong selection)
kansl1_calls <- read_tsv("clues/inv_dup_anc/rs80028338-calls.tsv", col_types = cols())

# load the filtered sample metadata
meta <- read_tsv("data/Ancestral_paths_new/ancestral_paths_merged_filtered_age.sampleInfo.tsv", col_types = cols())

inv_dup_calls <- dup_calls %>%

  # only retain West Eurasian samples from the main selection analysis
  inner_join(meta, by = c("Sample" = "sampleId")) %>%

  # join Andres' larger set of inversion calls
  left_join(inv_calls, by = c("Sample" = "sample")) %>%

  # drop any disagreements (ignoring NA)
  filter(Inversion == andres_call | is.na(Inversion) | is.na(andres_call)) %>%

  # get the union of the two sets of calls
  mutate(inv_call = coalesce(Inversion, andres_call)) %>%

  # drop any missing calls
  drop_na(inv_call, H1D, H2D) %>%

  # sort by age
  arrange(age) %>%

  # only keep the necessary columns
  select(Sample, Population, age, inv_call, H1D, H2D) %>%

  # convert to a more usable format
  mutate(call = case_when(
    inv_call == "H1H1" & H1D == 0 & H2D == 0 ~ "H1_nodup,H1_nodup",
    inv_call == "H1H1" & H1D == 1 & H2D == 0 ~ "H1_nodup,H1_dup",
    inv_call == "H1H1" & H1D == 2 & H2D == 0 ~ "H1_dup,H1_dup",

    inv_call == "H1H2" & H1D == 0 & H2D == 0 ~ "H1_nodup,H2_nodup",
    inv_call == "H1H2" & H1D == 0 & H2D == 1 ~ "H1_nodup,H2_dup",
    inv_call == "H1H2" & H1D == 1 & H2D == 0 ~ "H1_dup,H2_nodup",
    inv_call == "H1H2" & H1D == 1 & H2D == 1 ~ "H1_dup,H2_dup",

    inv_call == "H2H2" & H1D == 0 & H2D == 0 ~ "H2_nodup,H2_nodup",
    inv_call == "H2H2" & H1D == 0 & H2D == 1 ~ "H2_nodup,H2_dup",
    inv_call == "H2H2" & H1D == 0 & H2D == 2 ~ "H2_dup,H2_dup",
  )) %>%

  # filter(is.na(call)) %>%
  # write_tsv("inv17/Ancient_17q_Inversion_Duplication_Pop_Age-outliers.tsv")

  # drop any unusual allelic states
  drop_na(call)

# split the diploid calls into haploid calls, for the conditional modelling
haploid_calls <- inv_dup_calls %>%

  # split the comma delimited calls into separate lines
  separate_rows(call, sep = ",") %>%

  # split the underscore delimited calls into separate columns
  separate(call, into = c("inv", "dup"), sep = "_") %>%

  # only keep the necessary columns
  select(Sample, Population, age, inv, dup)

# convert the H1 vs H2 inversion calls into haploid CLUES input format
haploid_calls %>%
  mutate(call = case_when(
    inv == "H1" ~ "0.000000 -inf",
    inv == "H2" ~ "-inf 0.000000"
  )) %>%
  mutate(clues = sprintf("%.6f %s", age / 28, call)) %>%
  select(clues) %>%
  write_delim("clues/inv_dup/chr17-H1_vs_H2-clues.ancient", col_names = F, quote_escape = F, delim = "")

# convert the H1 ancestral vs duplication calls into haploid CLUES input format
haploid_calls %>%
  filter(inv == "H1") %>%
  mutate(call = case_when(
    dup == "nodup" ~ "0.000000 -inf",
    dup == "dup" ~ "-inf 0.000000"
  )) %>%
  mutate(clues = sprintf("%.6f %s", age / 28, call)) %>%
  select(clues) %>%
  write_delim("clues/inv_dup/chr17-H1_nodup_vs_dup-clues.ancient", col_names = F, quote_escape = F, delim = "")

# convert the H2 ancestral vs duplication calls into haploid CLUES input format
haploid_calls %>%
  filter(inv == "H2") %>%
  mutate(call = case_when(
    dup == "nodup" ~ "0.000000 -inf",
    dup == "dup" ~ "-inf 0.000000"
  )) %>%
  mutate(clues = sprintf("%.6f %s", age / 28, call)) %>%
  select(clues) %>%
  write_delim("clues/inv_dup/chr17-H2_nodup_vs_dup-clues.ancient", col_names = F, quote_escape = F, delim = "")


haplotypes <- c("H1_nodup", "H1_dup", "H2_nodup", "H2_dup")

for (hap in haplotypes) {
  # convert the inversion/duplication and KANSL1 calls into diploid CLUES input format
  inv_dup_calls %>%

    # merge the KANSL1 genotype calls
    inner_join(kansl1_calls, by = c("Sample" = "sample")) %>%

    # now split the inversion / duplication calls -- because we don't know the true phase of these :(
    separate_rows(call, sep = ",", convert = TRUE) %>%

    # split the GP calls into their three states
    separate(GP, sep = ",", into = c("HOM_REF", "HET", "HOM_ALT"), convert = TRUE, remove = FALSE) %>%

    # filter for the current inversion haplotype
    filter(call == hap) %>%

    # convert to CLUES format
    mutate(clues = sprintf("%.6f %.6f %.6f %.6f", age / 28, log10(HOM_REF), log10(HET), log10(HOM_ALT))) %>%
    select(clues) %>%
    write_delim(paste0("clues/inv_dup/chr17-", hap, "-KANSL1-clues.ancient"), col_names = F, quote_escape = F, delim = "")
}

# convert the KANSL1 genotype calls into haploid calls, so we can partition by ancestry
kansl1_haps <- inv_dup_calls %>%

  # merge the KANSL1 genotype calls
  inner_join(kansl1_calls, by = c("Sample" = "sample")) %>%

  # split the genotypes and ancestry paths into separate rows
  separate_rows(GT, AP, convert = TRUE) %>%

  # and split the basal ancestry paths
  mutate(AP = case_when(
    AP == 5 ~ "3|4",
    AP == 6 ~ "1|2",
    TRUE ~ as.character(AP)
  )) %>%
  separate_rows(AP, convert = TRUE) %>%

  # now split the inversion / duplication calls -- because we don't know the true phase of these :(
  separate_rows(call, sep = ",", convert = TRUE)

ancesties <- c("ANA" = 1, "CHG" = 2, "WHG" = 3, "EHG" = 4)

for (hap in haplotypes) {
  for (ap in names(ancesties)) {
    # CLUES input format
    kansl1_haps %>%
      filter(call == hap & AP == ancesties[ap]) %>%
      mutate(call = case_when(
        GT == 0 ~ "0.000000 -inf",
        GT == 1 ~ "-inf 0.000000"
      )) %>%
      mutate(clues = sprintf("%.6f %s", age / 28, call)) %>%
      select(clues) %>%
      write_delim(paste0("clues/inv_dup_anc/chr17-", hap, "-", ap, "-clues.ancient"), col_names = F, quote_escape = F, delim = "")
  }
}

inv_dup_calls %>%

  # merge the KANSL1 genotype calls
  inner_join(kansl1_calls, by = c("Sample" = "sample")) %>%

  # split the genotypes and ancestry paths into separate rows
  separate_rows(GT, AP, convert = TRUE) %>%
  group_by(inv_call) %>%
  summarise(DAC = sum(GT), count = n()) %>%
  mutate(DAF = DAC / count) %>%
  arrange(desc(DAF))

inv_dup_calls %>%

  # merge the KANSL1 genotype calls
  inner_join(kansl1_calls, by = c("Sample" = "sample")) %>%

  # split the genotypes and ancestry paths into separate rows
  separate_rows(GT, AP, convert = TRUE) %>%
  group_by(GT) %>%
  summarise(H1D_sum = sum(H1D), H2D_sum = sum(H2D), count = n())
