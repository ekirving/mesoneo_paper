#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(readr))
quiet(library(dplyr))
quiet(library(tidyr))

# load the filtered sample metadata
meta <- read_tsv("data/ancestral_paths_v3/ancestral_paths_v3.sampleInfo.tsv", col_types = cols()) %>%
    mutate(gens = age / 28)

# load the APOE SNPs
snps <- read_tsv("APOE.snps", col_names=c("ID", "REF", "ALT", "AA", "SAMPLE", "GT"), col_types = cols(.default = "c"))

data <- snps %>%
    separate(GT, into=c("h1", "h2"), convert=TRUE) %>%
    mutate(g1=ifelse(h1 == 0, REF, ALT), g2=ifelse(h2 == 0, REF, ALT)) %>%
    select(SAMPLE, ID, g1, g2) %>%
    pivot_longer(cols=c(g1, g2), names_to = "haplotype") %>%
    pivot_wider(id_cols = c("SAMPLE", "haplotype"), names_from="ID") %>%
    inner_join(meta %>% select(sampleId, gens), by = c("SAMPLE" = "sampleId")) %>%
    arrange(gens) %>%
    mutate(apoe = case_when(
        rs429358 == "T" & rs7412 == "T" ~ "APOE2",
        rs429358 == "T" & rs7412 == "C" ~ "APOE3",
        rs429358 == "C" & rs7412 == "C" ~ "APOE4",
        TRUE ~ "???"
    ))

suppressWarnings(dir.create(file.path("apoe")))

for (isoform in c("APOE2", "APOE3", "APOE4")) {
    
    # make the CLUES input file
    data %>%
        filter(gens != 0) %>%
        mutate(call = ifelse(apoe == isoform, "-inf 0.000000", "0.000000 -inf ")) %>%
        mutate(clues = sprintf("%.6f %s", gens, call)) %>%
        select(clues) %>%
        write_delim(paste0("apoe/", isoform, ".ancient"), col_names = F,  delim = "")
    
    # determine the modern frequency
    mod <- data %>% filter(gens == 0)
    freq <- sum(mod$apoe == isoform) / nrow(mod)
    cat(freq, file=paste0("apoe/", isoform, ".freq"), sep="\n")
}
