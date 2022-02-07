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
quiet(library(data.table))

# get the command line arguments
p <- arg_parser("Prepare the CLUES data for peak detection with harvester")
p <- add_argument(p, "--data", help = "CLUES data")
p <- add_argument(p, "--pairs", help = "GWAS and Control SNP pairings")
p <- add_argument(p, "--unmapped", help = "Unmappable SNPs from Relate")
p <- add_argument(p, "--flipped", help = "Flipped SNPs from Relate")
p <- add_argument(p, "--inlimit", help = "Filter for harvester input", default = 0.001)
p <- add_argument(p, "--output", help = "Output file")

argv <- parse_args(p)

# load all the data
data <- fread(argv$data, header = T, sep = "\t")
pairs <- fread(argv$pairs, header = T, sep = "\t")
unmapped <- fread(argv$unmapped, header = T, sep = "\t")
flipped <- fread(argv$flipped, header = T, sep = "\t")

# turn off the new summarise warnings
options(dplyr.summarise.inform = FALSE)

data <- data %>%

  # sort and add a bed format column
  arrange(chrom, start) %>%
  mutate(pos = paste0(chrom, ":", start)) %>%

  # drop any modern SNPs that were flipped or unmapped
  filter(!(mode == "Modern" & (pos %in% unmapped$pos | pos %in% flipped$pos))) %>%

  # only keep GWAS SNPs
  filter(rsid %in% pairs$gwas) %>%

  # harvester expects dense sampling, so compact the SNP positions
  mutate(adjusted = round(start / 10))

# harvester crashes with a floating point exception if there are chroms with no retained SNPs
bad <- data %>%
  group_by(chrom) %>%
  summarise(min = min(p.value)) %>%
  filter(min >= argv$inlimit) %>%
  pull(chrom)

# drop the bad chroms
data <- data %>%
  filter(!chrom %in% bad) %>%
  select(chrom, adjusted, p.value)

# save the data
fwrite(data, argv$output, sep = "\t", quote = FALSE)
