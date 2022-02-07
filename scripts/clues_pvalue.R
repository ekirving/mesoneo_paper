#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(data.table))

# get the command line arguments
p <- arg_parser("Calculare p-values for CLUES data")
p <- add_argument(p, "--data", help = "CLUES data")
p <- add_argument(p, "--output", help = "Output file")

argv <- parse_args(p)

# load all the data
data <- fread(argv$data, header = T, sep = "\t")

# calculate p-values from the log-likelihood ratio
data$p.value <- formatC(pchisq(data$logLR, df = 1, lower.tail = FALSE), format = "e", digits = 4)

# make sure the table is in order
data <- data[order(chrom, start)]

# save the data.table
fwrite(data, argv$output, sep = "\t", quote = FALSE)
