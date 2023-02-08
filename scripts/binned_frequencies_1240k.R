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
quiet(library(jsonlite))

# get the command line arguments
p <- arg_parser("Plot the binned allele frequencies for a SNP")
p <- add_argument(p, "--samples", help = "Sample metadata", default = "data/1240k/Le_et_al_2022_Table_S1.txt")
p <- add_argument(p, "--geno", help = "Pseudohaploid genotypes", default = "binned/1240k-lct.tsv")
p <- add_argument(p, "--output", help = "PNG file to output", default = "binned/lct-binned-calls-1240k.png")

argv <- parse_args(p)

# load the sample metadata
samples <- read_tsv(argv$samples, guess_max = 100000, col_types = cols()) %>%
    select(
        # tidy up the columns names
        sampleId=`Sample_ID`,
        ageAverage=`Date mean in BP in years before 1950 CE [OxCal mu for a direct radiocarbon date, and average of range for a contextual date]`
    )

# load the genotype calls
geno <- read_tsv(argv$geno, col_types = cols(.default = "c")) %>%
    # join the sample metadata
    inner_join(samples, by="sampleId") %>%
    # pivot SNPs from cols to rows
    pivot_longer(cols=-c("fid", "sampleId", "father", "mother", "sex", "phenotype", "ageAverage"), names_to = "rsid", values_to = "call") %>%
    # drop samples with no call
    filter(call != "0 0") %>%
    # split calls
    separate(call, into = c("h1", "h2")) %>%
    pivot_longer(cols=c("h1", "h2"), names_to = "haplotype", values_to = "genotype")

# get the ancestral allele for all SNPs
anc <- data.frame()
for (rsid in unique(geno$rsid)) {
    meta <- fromJSON(paste0("variants/metadata/GRCh37/", rsid, ".json"))
    anc <- bind_rows(anc, data.frame("rsid"=c(rsid), "ancestral"=c(meta$ancestral)))
}

# polarize the genotypes
geno <- geno %>%
    inner_join(anc) %>%
    mutate(derived=as.numeric(genotype != ancestral))

# determine the sample size of the ancients
sample_size <- geno %>% select(sampleId) %>% unique() %>% nrow()

BIN_SIZE <- 1000

# group by age
binned <- geno %>%
    # add a new age bin
    mutate(bin=-round(ageAverage/BIN_SIZE)*BIN_SIZE) %>%
    group_by(rsid, bin) %>%
    summarise(sum=sum(derived), count=n(), daf=sum/count)

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

    guides(size=guide_legend(title="Genotype count")) +

    scale_x_continuous(breaks = xbreaks, labels = xlabels) +
    scale_size_continuous(limits  = c(size_min, size_max), breaks = size_breaks) +
    labs(title = paste0("v52.2_1240K_public (n=", sample_size," ancient samples)")) +

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

ggsave(argv$output, width = 5.5, height = 4)

