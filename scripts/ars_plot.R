#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2022, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
    suppressMessages(suppressWarnings(x))
}
quiet(library(argparser)) # v0.6
quiet(library(tidyverse)) # v1.3.1

# load the ARS data 
data <- read_csv("~/Downloads/ars-plot/processed_PRS_calculations_1000bootstrapped_v3_split_192001-216000_phenotypes_alba_5e-8_lowestpval.csv", col_types = cols()) %>%
    mutate(phenocode=str_replace(str_replace(phenotype, "/willerslev/datasets/UKBiobank/NealeV2/", ""), ".gwas.imputed_v3.both_sexes.tsv.gz", "")) %>%
    select(phenocode, ancestry, mean, std)

# load the manifest
maifest <- read_csv("~/Downloads/ars-plot/Manifest_201807.csv", col_types = cols(), na = c("", "NA", "-")) %>%
    filter(`Sex`=="both_sexes") %>%
    select(phenocode=`Phenotype Code`, name=`Phenotype Description`) %>%
    mutate(name = str_replace(name, "Diagnoses - main ICD10: ", ""))

# join the trait names to the codes
data <- data %>%
    inner_join(maifest, by = 'phenocode', relationship = "many-to-one") %>%
    drop_na()

top_ancestry <- "Farmer"

# set the sort order of the marginal traits based on their R-score in the ancestry with the largest score
data$name <- factor(data$name, levels = data %>% filter(ancestry == top_ancestry) %>% arrange(mean) %>% pull(name) %>% unique())

data$ancestry <- factor(data$ancestry, levels = c("WHG", "EHG", "CHG", "Farmer", "Yamnaya"), labels = c("WHG", "EHG", "CHG", "Farmer"="Neolithic Farmer", "Yamnaya"="Steppe"))

ancestry_colors <- c(
    "WHG" = "#fc8d62",
    "EHG" = "#8da0cb",
    "CHG" = "#e78ac3",
    "Neolithic Farmer" = "#a6d854",
    "Steppe" = "#b15928"
)
                        
data %>%
    ggplot(aes(x = name, y = mean, color = ancestry)) +
    # add a solid line at the zero mark
    geom_hline(yintercept = 0, color = "darkgrey") +
    geom_errorbar(aes(ymin=mean-(std*1.96), ymax=mean+(std*1.96)), width=.3)+
    ylab("Ancestral Risk Score (z-score)") +
    
    # set the colour scales
    scale_color_manual(values = ancestry_colors) +
    
    # basic styling
    theme_bw() +
    theme(
        # legend.position = "none",
        legend.position="top",
        legend.title=element_blank(),
        legend.text = element_text(size = 7),
        legend.margin=margin(c(0,0,0,0)),
        # legend.position="bottom",
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "white", color = "white"),
        plot.margin = margin(0.1, 0.1, 0, 0.5, "cm"),
        panel.spacing = unit(0.5, "lines"),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1)
    )

# save the plot
ggsave("figs_hires/Figure_6.pdf", height = 10, width = 9, units = "cm", scale=1.4)

