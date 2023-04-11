#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(dplyr))
quiet(library(stringr))
quiet(library(tidyr))
quiet(library(readr))

# get the command line arguments
p <- arg_parser("Detect selective sweeps in the CLUES data")
p <- add_argument(p, "--data", help = "CLUES data", default = "clues/ancestral_paths_v3-all-ancient-ALL-filtered-clues_report.tsv")
p <- add_argument(p, "--pairs", help = "GWAS and neutral SNP pairings", default = "variants/ancestral_paths_v3-all-pairs.tsv")
p <- add_argument(p, "--unmapped", help = "Unmappable SNPs from Relate", default = "relate/1000G_phase3-FIN_GBR_TSI-popsize-allsnps_unmapped.tsv.gz")
p <- add_argument(p, "--flipped", help = "Flipped SNPs from Relate", default = "relate/1000G_phase3-FIN_GBR_TSI-popsize-allsnps_flipped.tsv.gz")
p <- add_argument(p, "--1240k", help = "Data is from the 1240k SNP array", flag = TRUE)
p <- add_argument(p, "--output", help = "List of peaks", default = "clues/ancestral_paths_v3-all-ancient-ALL-filtered-sweeps.tsv")

argv <- parse_args(p)

# load all the data
data <- read_tsv(argv$data, col_types = cols())

if (!argv$`1240k`) {
    # minimum number of SNPs to call a sweep
    MIN_SNPS <- 6
    
    pairs <- read_tsv(argv$pairs, col_types = cols())
    unmapped <- read_tsv(argv$unmapped, col_types = cols(.default = "c"))
    flipped <- read_tsv(argv$flipped, col_types = cols(.default = "c"))
    
    data <- data %>%
        # determine if SNPs are GWAS or controls
        mutate(type = ifelse(rsid %in% pairs$gwas, "gwas", "control")) %>%
    
        # drop any modern SNPs that were flipped or unmapped
        mutate(pos = paste0(chrom, ":", start)) %>%
        filter(!(mode == "modern" & (pos %in% unmapped$pos | pos %in% flipped$pos)))
} else {
    # use a smaller threshold for the 1240k data
    MIN_SNPS <- 3
    
    data <- data %>% 
        mutate(chrom=CHR, start=POS, end=POS, p.value=corrected.p, ancestry=NA, type="1240k")
}

# use a genome-wide significance threshold
p.genomewide <- 5e-8

cluster_chrom <- function(data, min_snps) {

    # only cluster significant SNPs
    data <- data %>%
        filter(p.value <= p.genomewide)

    # skip chromosomes with an insufficient number of SNPs
    if (nrow(data) < 2) {
        return(NULL)
    }

    # calculate the distance between all significant SNPs
    dist_mat <- dist(data %>% select(start))

    # hierarchically cluster the SNPs, using the 'single linkage' method
    hclust_avg <- hclust(dist_mat, method = "single")

    # split the clusters on a max height of 1 Mb
    data$cluster <- cutree(hclust_avg, h = 1e6)

    data %>%
        group_by(cluster) %>%
        mutate(
            # get the co-ordinates of the cluster
            start = min(start),
            end = max(end),
            # label the clusters by their first and last SNP
            label = sprintf(
                "chr%s:%.1f-%.1f Mb",
                chrom,
                start / 1e6,
                end / 1e6
            ),
            # count the SNPs in each cluster
            num_snps = n()
        ) %>%
        slice_min(p.value, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        select(type, ancestry, label, chrom, start, end, num_snps, top_snp=rsid, min_pval=p.value) %>%
        # drop clusters below the minimum threshold
        filter(num_snps >= min_snps)
}

# cluster the GWAS and Control SNPs separately for each chromosome
sweeps <- bind_rows(
    lapply(unique(data$type), function(snp_type) {
        group <- data %>% filter(type == snp_type)
        bind_rows(
            lapply(seq(22), function(chr) {
                cluster_chrom(group %>% filter(chrom == chr), MIN_SNPS)
            })
        )
    })
)


if (nrow(sweeps) == 0) {
  # handle edge case of insufficient SNPs to cluster
  sweeps <- tibble(type=character(), ancestry=character(), label=character(), chrom=character(), start=numeric(), end=numeric(), num_snps=numeric(), top_snp=character(), min_pval=numeric())
}
   
write_tsv(sweeps, argv$output)