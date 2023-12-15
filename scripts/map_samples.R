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
quiet(library(sf))
quiet(library(rnaturalearth))
quiet(library(rnaturalearthdata))
quiet(library(rgeos))
quiet(library(ggthemes))
quiet(library(viridis))
quiet(library(tidyquant))
quiet(library(ggdist))
quiet(library(gghalves))
quiet(library(ggpubr))

# get the command line arguments
p <- arg_parser("Plot the map of sampling locations")
p <- add_argument(p, "--samples", help = "Sample metadata", default = "data/ancestral_paths_v3/ancestral_paths_v3.sampleInfo.tsv")
p <- add_argument(p, "--output", help = "PDF file to output", default = "figs_hires/Figure_1.pdf")

argv <- parse_args(p)

region_levels <- c("CentralWesternAsia", "NorthernEurope", "SouthernEurope", "CentralEasternEurope", "WesternEurope")
region_labels <- c("Central / Western Asia", "Northern Europe", "Southern Europe", "Central / Eastern Europe", "Western Europe")

# load the sample metadata
samples <- read_tsv(argv$samples, col_types = cols(), guess_max=5000) %>%
    # drop the modern samples
    filter(groupAge != "Modern") %>%
    # fix the issue with longitude
    mutate(longitude=as.numeric(trimws(longitude, whitespace = "[\\h\\v]"))) %>%
    # merge Central and Western Asia
    mutate(region = ifelse(region %in% c("CentralAsia", "WesternAsia"), "CentralWesternAsia", region)) %>%
    # convert to factors so we can set the display order
    mutate(region = factor(region, levels=region_levels, labels=region_labels)) %>%
    # only keep the necessary columns
    select(sampleId, popId, region, ageAverage, latitude, longitude) %>%
    # sort by age so the oldest are printed on top
    arrange(ageAverage)

# convert to an `sf` object
samples_sf <- st_as_sf(samples, coords=c('longitude', 'latitude'), crs=4326)

# fetch the map
world <- ne_countries(scale = "medium", returnclass = "sf")
europe <- ne_countries(continent = "europe", returnclass = "sf")

plt_map <- ggplot() +
    # use a world map
    geom_sf(data = world) +
    # display the samples/ages
    geom_sf(data = samples_sf, mapping = aes(color=ageAverage), position="dodge") +
    # crop the map to only show West Eurasia
    coord_sf(xlim = c(-25, 77), ylim = c(33, 71.5), expand = FALSE) +
    # mask the edges of the map
    annotate("rect", xmin = -25, xmax = -0, ymin = 68, ymax = 100, fill = "white") +
    scale_color_viridis_c(limits = c(0, 15000), option="magma") +
    labs(colour="Sample Age") +
    theme_map() +
    theme(plot.background = element_rect(fill = "white", color="white"))

# -----------------------------------------------------------------------------------------------

# constrain the extent of the plotting
xbreaks <- seq(0, 15000, 1000)
xlabels <- (xbreaks / 1000)
        
plt_ages <- samples %>%
    ggplot(mapping=aes(x = region, y = ageAverage, fill = factor(region))) +
    # add half-violin from {ggdist} package
    stat_halfeye(
        # adjust bandwidth
        adjust = 0.5,
        # move to the right
        justification = -0.2,
        # remove the slub interval
        .width = 0,
        point_colour = NA
    ) +
        geom_boxplot(
            width = 0.12,
            # removing outliers
            outlier.color = NA,
            alpha = 0.5
        ) +
        stat_dots(
            # ploting on left side
            side = "left",
            # adjusting position
            justification = 1.1,
            # adjust grouping (binning) of observations
            binwidth = 50
            # color="black"
        ) +
    scale_fill_tq() +
    theme_tq() +
    theme(
        axis.title.y=element_blank(),
        axis.text.y = element_text(vjust=-3, hjust=0, size=9, margin = margin(0, -105, 0, 0)),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none"
    ) +
    scale_y_reverse(breaks = xbreaks, labels = xlabels) +
    ylab("kyr BP") +
    coord_flip()

# -----------------------------------------------------------------------------------------------

# arrange the columns into a column
main_fig <- ggarrange(
    plotlist = list(plt_map, plt_ages),
    labels = c("a", "b"),
    label.x = c(0.03, -0.03),
    ncol = 2,
    nrow = 1,
    widths=c(4,3)
)

ggsave(filename = argv$output, plot = main_fig, width=11, height=3.8, limitsize = FALSE)