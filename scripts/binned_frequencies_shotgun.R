#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(tidyverse))
quiet(library(directlabels))

# get the command line arguments
p <- arg_parser("Plot the binned allele frequencies for a SNP")
p <- add_argument(p, "--samples", help = "Sample metadata", default = "data/Ancestral_paths_new/ancestral_paths_merged_filtered.sampleInfo.tsv")
p <- add_argument(p, "--imputed", help = "Imputed genotypes", default = "binned/ancestral_paths_new-lct-genotypes.tsv")
p <- add_argument(p, "--likelihoods", help = "Genotype likelihoods", default = "binned/neo_likelihoods-lct-genotypes.tsv")
p <- add_argument(p, "--output", help = "PNG file to output", default = "binned/lct-binned-calls-shotgun.png")

argv <- parse_args(p)

# load the sample metadata
samples <- read_tsv(argv$samples, col_types = cols()) %>%
    select(sampleId, popId, groupAge, ageAverage)

# load the genotype calls
impute <- read_tsv(argv$imputed, col_types = cols(.default = "c"), col_names = c("rsid", "ref", "alt", "anc", "sample", "geno"))
gl <- read_tsv(argv$likelihoods, col_types = cols(.default = "c"), col_names = c("rsid", "ref", "alt", "anc", "sample", "geno"))

# merge the two callsets
geno <- bind_rows(
    gl %>% mutate(type="gl"),
    impute %>% mutate(type="impute")
) %>%
    separate(col="rsid", sep = ";", into = c("rsid", "pos"), fill="right") %>%
    # drop missing genotypes
    filter(!geno %in% c(".", "./.")) %>%
    # split "0/0" and "0|0" notation into separate rows
    separate(col="geno", into=c("hap1_geno", "hap2_geno"), sep="[|/]") %>%
    pivot_longer(cols=starts_with("hap"), names_to = c("haplotype", ".value"), names_sep="_") %>%
    # convert numeric calls into alleles
    mutate(
        call=ifelse(geno==0, ref, alt),
        derived=ifelse(call==anc, 0, 1)
    ) %>%
    # join the sample metadata
    inner_join(samples, by=c("sample"="sampleId"))

# determine the sample size of the ancients
sample_size <- geno %>% filter(groupAge == "Ancient") %>% select(type, sample) %>% unique() %>% group_by(type) %>% tally()
sample_size_gl <- sample_size %>% filter(type == "gl") %>% pull(n)
sample_size_impute <- sample_size %>% filter(type == "impute") %>% pull(n)

BIN_SIZE <- 1000

# group by age
binned <- geno %>%
    # fill in missing ages for the modern samples
    mutate(ageAverage = ifelse(groupAge=="Modern", 0, ageAverage)) %>%
    # add a new age bin
    mutate(bin=-round(ageAverage/BIN_SIZE)*BIN_SIZE) %>%
    group_by(rsid, type, bin) %>%
    summarise(sum=sum(derived), count=n(), daf=sum/count) %>%
    mutate(type=ifelse(type=="gl",
                       paste0("Non-imputed calls (n=", sample_size_gl," ancient samples)"),
                       paste0("Imputed calls (n=", sample_size_impute," ancient samples)")))

# constrain the extent of the plotting
xmin <- min(binned$bin)
xmax <- max(binned$bin)
xbreaks <- seq(xmin, xmax, BIN_SIZE)
xlabels <- round(xbreaks / BIN_SIZE)

size_min <- floor(min(binned$count) / 100) * 100
size_max <- ceiling(max(binned$count) / 100) * 100
size_breaks <- 2^(1:ceiling(log2(size_max)))

ggplot(binned, aes(x=bin, y=daf, color=rsid)) +

    geom_smooth(method = "loess", formula = "y ~ x", se = FALSE) +
    geom_point(aes(size=count), alpha=.8, position=position_dodge(width = 300)) +
    facet_wrap(~type) +

    guides(size=guide_legend(title="Genotype count")) +

    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    scale_size_continuous(limits  = c(size_min, size_max), breaks = size_breaks) +

    ylab("DAF") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
        # legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        panel.spacing = unit(1, "lines")
    )

ggsave(argv$output, width = 11, height = 4)

