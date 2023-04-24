#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(readr))
quiet(library(bedr))
quiet(library(dplyr))
quiet(library(tidyr))

# # get the command line arguments
p <- arg_parser("Plot the main text CLUES figure")
p <- add_argument(p, "--peaks-all", help = "CLUES selection peaks from the pan-ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-ALL-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-whg", help = "CLUES selection peaks from the WHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-WHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-ehg", help = "CLUES selection peaks from the EHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-EHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-chg", help = "CLUES selection peaks from the CHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-CHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-ana", help = "CLUES selection peaks from the ANA ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-ANA-filtered-sweeps.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "clues/ancestral_paths_v3-all-ancient-filtered-sweeps-merged.tsv")

argv <- parse_args(p)

# load some helper functions
source("scripts/clues_utils.R")

# load all the CLUES peaks (across all ancestries)
ancient_peaks <- bind_rows(
    read_tsv(argv$peaks_all, col_types = cols(.default="n", type="c", ancestry="c", label="c", top_snp="c")),
    read_tsv(argv$peaks_whg, col_types = cols(.default="n", type="c", ancestry="c", label="c", top_snp="c")),
    read_tsv(argv$peaks_ehg, col_types = cols(.default="n", type="c", ancestry="c", label="c", top_snp="c")),
    read_tsv(argv$peaks_chg, col_types = cols(.default="n", type="c", ancestry="c", label="c", top_snp="c")),
    read_tsv(argv$peaks_ana, col_types = cols(.default="n", type="c", ancestry="c", label="c", top_snp="c")),
)

if (nrow(ancient_peaks) > 0) {
    # merge overlapping sweep regions
    gwas <- sweep_to_bed(ancient_peaks %>% filter(type == "gwas"))
    control <- sweep_to_bed(ancient_peaks %>% filter(type == "control"))
    
    # convert from bed format back into a table
    sweeps <- bind_rows(
        tibble(type=c("gwas"), bed=gwas),
        tibble(type=c("control"), bed=control)
    ) %>% 
        separate(col = bed, into = c("chrom", "start", "end"), sep = "[:-]", convert = TRUE) %>%
        mutate(
            # strip the chr prefix
            chrom=as.integer(sub("chr", "", chrom)),
            # label the clusters by their first and last SNP
            label = sprintf(
                "chr%s:%.1f-%.1f Mb",
                chrom,
                start / 1e6,
                end / 1e6
            ),
            .before=chrom
        ) %>%
        arrange(chrom, start)
} else {
    # handle edge case of no sweeps to merge
    sweeps <- tibble(type=character(), label=character(), chrom=character(), start=numeric(), end=numeric())
}

write_tsv(sweeps, argv$output)
