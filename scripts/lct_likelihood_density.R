#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(tidyverse))

# load the sample metadata
samples <- read_tsv("data/ancestral_paths_v3/ancestral_paths_merged_filtered.sampleInfo.tsv", col_types = cols()) %>%
    select(sampleId, ageAverage)

# load the GLs
data <- read_tsv("lct-sweep/lct-likelihoods.tsv", guess_max = 100000, col_types = cols(),
                 col_names = c("rsid", "ref", "alt", "ancestral", "sampleId", "genotype", "likelihood")) %>%
    # drop samples with no call
    filter(likelihood != "." & genotype != "./.") %>%
    # drop unneeded columns
    select(rsid, sampleId, likelihood) %>%
    # join the sample age
    inner_join(samples, by="sampleId") %>%
    # split the GLs into separate columns
    separate(likelihood, into = c("hom_ref", "het", "hom_alt"), convert = TRUE) %>%
    # pivot longer
    pivot_longer(cols=c("hom_ref", "het", "hom_alt"), names_to="genotype", values_to="pl") %>%
    # convert from Phred-scaled likelihoods back to conditional probabilities (i.e, `P(Genotype | Data)`)
    # see https://gatk.broadinstitute.org/hc/en-us/articles/360035890451-Calculation-of-PL-and-GQ-by-HaplotypeCaller-and-GenotypeGVCFs
    mutate(cond_prob=10 ** (-pl / 10)) %>%
    # de-normalize
    group_by(rsid, sampleId) %>%
    mutate(cond_prob = cond_prob/sum(cond_prob))

# set facet ordering
data$genotype <- factor(data$genotype, levels = c("hom_ref", "het", "hom_alt"), labels = c("hom-ref", "het", "hom-alt"))

break_size <- 60
xmin <- min(data$pl)
xmax <- ceiling(max(data$pl) / break_size) * break_size
xbreaks <- seq(xmin, xmax, break_size)

# plot the Phred-scaled likelihoods
ggplot(data, aes(x=pl, group=genotype, fill=genotype)) +
    geom_histogram() +
    facet_grid(rsid ~ genotype) +
    scale_x_continuous(breaks = xbreaks) +

    ylab("Genotype count") +
    xlab("Phred scaled likelihoods") +

    # basic styling
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines")
    )

ggsave("lct-sweep/lct-likelihoods.png", width = 6, height = 4)

# plot the conditional probabilities
ggplot(data, aes(x=cond_prob, group=genotype, fill=genotype)) +
    geom_histogram() +
    facet_grid(rsid ~ genotype) +

    ylab("Genotype count") +
    xlab("P(Genotype | Data)") +

    # basic styling
    theme_minimal() +
    theme(
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines")
    )

ggsave("lct-sweep/lct-cond_prob.png", width = 6, height = 4)
