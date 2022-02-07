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
quiet(library(ggplot2))
quiet(library(RcppCNPy))
quiet(library(RcppRoll))
quiet(library(tibble))
quiet(library(tidyr))

# get the command line arguments
p <- arg_parser("Plot CLUES trajectory with ggplot")
p <- add_argument(p, "--prefix", help = "Prefix for the CLUES model data")
p <- add_argument(p, "--gen-time", default = 28, help = "Generation time")
p <- add_argument(p, "--max-age", default = 13665, help = "Upper bound for plotting")
p <- add_argument(p, "--output", help = "Output file")

argv <- parse_args(p)

# load the model data
epochs <- npyLoad(paste0(argv$prefix, ".epochs.npy"))
freqs <- npyLoad(paste0(argv$prefix, ".freqs.npy"))
logpost <- npyLoad(paste0(argv$prefix, ".post.npy"), dotranspose = F)

# constrain the extent of the plotting
xmin <- min(epochs)
xmax <- min(max(epochs), round(argv$max_age / argv$gen_time))
xbreaks <- seq(xmin, xmax + 1, round(1000 / argv$gen_time))
xlabels <- round(xbreaks * argv$gen_time / 1000)

df <- as_tibble(logpost) %>%

  # convert posterior densities from log-likelihoods
  exp() %>%

  # add the frequency labels
  add_column(freq = freqs, .before = 1) %>%

  # add the title heights (and a little padding to make sure there are no gaps)
  # tile heights are not equal because there is higher sampling desnity near 0 and 1
  add_column(height = diff(c(0, freqs)) + 1e-4, .before = 2) %>%

  # pivot the columns into long format
  pivot_longer(-c(freq, height), names_to = "epoch", values_to = "density") %>%

  # convert the column names into epochs (and switch the direction of time)
  mutate(epoch = -epochs[as.numeric(gsub("V", "", epoch))])

max_traj <- df %>%
  group_by(epoch) %>%
  top_n(1, density) %>%
  ungroup() %>%
  arrange(epoch) %>%
  mutate(freq = roll_mean(freq, 5, align = "left", fill = max(freq)))


plt <- df %>%
  # plot the heatmap
  ggplot(aes(x = epoch, y = freq, height = height, fill = density)) +
  geom_tile() +

  # plot the maximum posterior trajectory
  # geom_line(aes(x = epoch, y = freq), max_traj) +

  # set the axis breaks
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = c(0, 0), position = "right") +
  scale_x_continuous(limits = c(-xmax, xmin), breaks = -xbreaks, labels = xlabels, expand = c(0, 0)) +
  scale_fill_viridis_c(limits = c(0, 0.5)) +
  # scale_fill_gradient(low = "white", high = "black", limits = c(0, 0.5)) +

  labs(title = argv$prefix) +
  ylab("Derived Allele Frequency") +
  xlab("kyr BP") +

  # basic styling
  theme_minimal() +
  theme(
    # legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()
  )

png(file = argv$output, width = 16, height = 9, units = "in", res = 300)
plt
dev <- dev.off()
