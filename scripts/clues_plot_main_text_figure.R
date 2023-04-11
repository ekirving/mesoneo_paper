#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2021, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(argparser))
quiet(library(readr))
quiet(library(ggpubr))
quiet(library(gtools))

# # get the command line arguments
p <- arg_parser("Plot the main text CLUES figure")
p <- add_argument(p, "--data", help = "CLUES report", default = "clues/ancestral_paths_v3-all-filtered-clues_report.tsv")
p <- add_argument(p, "--pairs", help = "GWAS and neutral SNP pairings", default = "variants/ancestral_paths_v3-all-pairs.tsv")
p <- add_argument(p, "--unmapped", help = "List of modern SNPs unmappable by Relate", default = "relate/1000G_phase3-FIN_GBR_TSI-popsize-allsnps_unmapped.tsv.gz")
p <- add_argument(p, "--flipped", help = "List of modern SNPs flipped by Relate", default = "relate/1000G_phase3-FIN_GBR_TSI-popsize-allsnps_flipped.tsv.gz")
p <- add_argument(p, "--peaks-all", help = "Manhattan Harvester peaks from the pan-ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-ALL-filtered-harvester.tsv")
p <- add_argument(p, "--peaks-whg", help = "Manhattan Harvester peaks from the WHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-WHG-filtered-harvester.tsv")
p <- add_argument(p, "--peaks-ehg", help = "Manhattan Harvester peaks from the EHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-EHG-filtered-harvester.tsv")
p <- add_argument(p, "--peaks-chg", help = "Manhattan Harvester peaks from the CHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-CHG-filtered-harvester.tsv")
p <- add_argument(p, "--peaks-ana", help = "Manhattan Harvester peaks from the ANA ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-ANA-filtered-harvester.tsv")
p <- add_argument(p, "--mathieson", help = "Manhattan Harvester peaks from Mathieson et al. 2015", default = "mathieson/mathieson-harvester.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "figs/main-fig-original.png")

argv <- parse_args(p)

# load some helper functions
source("scripts/clues_utils.R")

ancestries <- c("ALL", "WHG", "EHG", "CHG", "ANA")

# load all the CLUES results
clues <- load_clues(argv$data, argv$pairs, argv$unmapped, argv$flipped)

# ----------------------------------------------------------------------------------------------------------------

# get the extra models for `rs1438307`
extra_snp <- read_tsv("rs1438307-report.tsv", col_types = cols()) %>%
    rename(rsid=rsid_x) %>%
    select(-c(derived, minor, rsid_y))

# calculate p-values from the log-likelihood ratio
extra_snp$p.value <- pchisq(2 * extra_snp$logLR, df = 1, lower.tail = FALSE)

# determine if SNPs are GWAS or controls
extra_snp$type <- "gwas"

# set facet ordering
extra_snp$mode <- factor(extra_snp$mode, levels = c("modern", "ancient"), labels = c("Modern", "Ancient"))
extra_snp$type <- factor(extra_snp$type, levels = c("gwas", "control"), labels = c("GWAS", "Control"))
extra_snp$ancestry <- factor(extra_snp$ancestry, levels = c("ALL", "WHG", "EHG", "CHG", "ANA"))

clues <- bind_rows(clues, extra_snp)

# ----------------------------------------------------------------------------------------------------------------

# Manhattan Harvester results for each of the CLUES analysis groups (i.e. modern, ancient, WHG, EHG, CHG and ANA)
ancient_all_peaks <- suppressWarnings(read_tsv(argv$peaks_all, col_types = cols()))
ancient_whg_peaks <- suppressWarnings(read_tsv(argv$peaks_whg, col_types = cols()))
ancient_ehg_peaks <- suppressWarnings(read_tsv(argv$peaks_ehg, col_types = cols()))
ancient_chg_peaks <- suppressWarnings(read_tsv(argv$peaks_chg, col_types = cols()))
ancient_ana_peaks <- suppressWarnings(read_tsv(argv$peaks_ana, col_types = cols()))

# Manhattan Harvester results for the reported p-values in Mathieson_et_al_2015
mathieson_peaks <- suppressWarnings(read_tsv(argv$mathieson, col_types = cols()))

# color brewer https://colorbrewer2.org/#type=qualitative&scheme=Set2&n=5
ancestry_colors <- c(
  "ALL" = "#66c2a5",
  "WHG" = "#fc8d62",
  "EHG" = "#8da0cb",
  "CHG" = "#e78ac3",
  "ANA" = "#a6d854"
)

# color brewer https://colorbrewer2.org/#type=qualitative&scheme=Paired&n=10
snp_colors <- c(
  "#1f78b4",
  "#33a02c",
  "#e31a1c",
  "#ff7f00",
  "#6a3d9a"
)

# get all the ancient ancestry models
ancestral <- filter(clues, mode == "Ancient")

# annotate ancestry specific peaks
peaks_all <- annotate_peaks(filter(ancestral, ancestry == "ALL"), ancient_all_peaks, mathieson_peaks)
peaks_whg <- annotate_peaks(filter(ancestral, ancestry == "WHG"), ancient_whg_peaks, mathieson_peaks)
peaks_ehg <- annotate_peaks(filter(ancestral, ancestry == "EHG"), ancient_ehg_peaks, mathieson_peaks)
peaks_chg <- annotate_peaks(filter(ancestral, ancestry == "CHG"), ancient_chg_peaks, mathieson_peaks)
peaks_ana <- annotate_peaks(filter(ancestral, ancestry == "ANA"), ancient_ana_peaks, mathieson_peaks)

# concatenate all the annotated peaks data back together
peaks_joined <- bind_rows(
  peaks_all$data,
  peaks_whg$data,
  peaks_ehg$data,
  peaks_chg$data,
  peaks_ana$data,
) %>% arrange(as.character(chrom), start)

# merge all the peak regions (across all ancestries)
all_regions <- havester_to_bed(
  bind_rows(
    ancient_all_peaks,
    ancient_whg_peaks,
    ancient_ehg_peaks,
    ancient_chg_peaks,
    ancient_ana_peaks,
  ) %>%
    # add 50kb buffer either side (for display purposes only; does not effect the reported intervals)
    separate(col = range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
    mutate(range = paste0(start - 5e5, "-", end + 5e5))
)

# join the merged regions
peaks_joined <- bedr.join.region(peaks_joined$bed, all_regions, verbose = FALSE) %>%
  filter(V4 != ".") %>%
  mutate(merged_peaks = paste0(V4, ":", as.numeric(V5) + 5e5, "-", as.numeric(V6) - 5e5)) %>%
  select(index, merged_peaks) %>%
  unique() %>%
  left_join(peaks_joined, ., by = c("bed" = "index"))


# for SNPs inside GWAS selection peaks, only retain the most significant marginal ancestry
ancestral_min <- peaks_joined %>%
  mutate(merged_peaks = ifelse(type == "Control", NA, merged_peaks)) %>%
  filter(type == "GWAS") %>%
  filter((ancestry == "ALL" & (type == "Control" | is.na(merged_peaks))) | (ancestry != "ALL" & (type == "GWAS" & !is.na(merged_peaks)))) %>%
  group_by(rsid) %>%
  slice_min(p.value, n = 1)

# get the gene names for each peak
gene_names <- peaks_joined %>%
  drop_na(merged_peaks) %>%
  separate_rows(genes, sep = " - ") %>%
  filter(genes != "" & genes != "NA") %>%
  group_by(chrom, merged_peaks, genes) %>%
  summarise(p.sum = sum(-log10(p.value)), .groups = "drop_last") %>%
  slice_max(p.sum) %>%
  separate(genes, c("genes", NA), sep = "[,;]", fill = "right") %>%
  mutate(genes = ifelse(grepl("chr6:.+", merged_peaks), "HLA", genes))

# use a genome-wide significance threshold
p.genomewide <- 5e-8

# count the number of SNPs in each group
num_gwas <- length(unique(filter(ancestral, type == "GWAS")$rsid))

# calculate the Bonferroni significance threshold
p.bonferroni <- 0.05 / num_gwas

# plot the ancestry stratified selection peaks
plt_man <- manhattan_plot(ancestral_min, c("type"), p.genomewide, p.bonferroni, show_lables = FALSE, show_legend = TRUE, show_strip = FALSE, composite = TRUE, colors = ancestry_colors, gene_names = gene_names) +
  theme(plot.title = element_text(hjust = 0.5))

# plt_man

# ------------------------------------------------------------------------------------------------------------------------------------------------------

# pick three regions for the zoomed-in columns
focal_regions <- gene_names %>%
  filter(genes %in% c("MCM6", "SLC45A2", "FADS2")) %>%
  pull(merged_peaks) %>%
  mixedsort()

# get the top SNP from each of the marginal ancestries
top_marginal_snps <- peaks_joined %>%
  drop_na(merged_peaks) %>%
  filter(p.value < p.bonferroni) %>%
  group_by(merged_peaks, ancestry) %>%
  slice_min(p.value, n = 1) %>%
  group_by(merged_peaks, rsid) %>%
  summarise(n = length(rsid), p.sum = sum(-log10(p.value))) %>%
  arrange(merged_peaks, desc(p.sum))

plt_columns <- lapply(focal_regions, function(focal_region) {

  # get all the SNPs in this region
  region_data <- filter(peaks_joined, merged_peaks == focal_region) %>%
      # remove rs12465802
      filter(rsid != "rs12465802")

  # get the list of top SNPs in the region
  region_snps <- filter(top_marginal_snps, merged_peaks == focal_region)$rsid

  # ----------------------------------------------------------------------------------------------------------------

  if (focal_region == "chr2:135300859-137564020") {
    # add the extra SNP
    region_snps <- c(region_snps, "rs1438307")
  }

  # ----------------------------------------------------------------------------------------------------------------

  # pair the colours with the ordered top SNPs
  label_colors <- setNames(snp_colors[1:length(region_snps)], region_snps)

  # make region name human readable
  region_name <- unlist(strsplit(focal_region, split = "[:-]"))
  region_name <- sprintf("%s:%.1f-%.1f Mb", region_name[1], as.integer(region_name[2]) / 1e6, as.integer(region_name[3]) / 1e6)

  # determine the max height, to keep the y-axis in syc
  max_height <- max(-log10(region_data$p.value))

  # make the zoomed in manhattan plots
  plt_col_zoom <- lapply(ancestries, function(focal_ancestry) {
    manhattan_plot(filter(region_data, ancestry == focal_ancestry), c("merged_peaks", "ancestry"), p.genomewide, p.bonferroni, top_snps = region_snps, show_strip = FALSE, wrap = 1, composite = FALSE, colors = ancestry_colors, label_colors = label_colors, region_name = region_name, max_height = max_height)
  })

  # make the per-ancestry CLUES plots
  plt_col_clues <- lapply(ancestries, function(focal_ancestry) {

    # sort the SNPs by their marginal p-value
    ordered_snps <- region_data %>%
      filter(rsid %in% region_snps & ancestry == focal_ancestry) %>%
      arrange(p.value) %>%
      pull(rsid)

    clues_plot("ancestral_paths_v3", "all", ordered_snps, "ancient", focal_ancestry, label_colors)
  })

  row_label <- filter(gene_names, merged_peaks == focal_region)$genes

  # arrange the plots into two columns
  ggarrange(
    ggarrange(plotlist = plt_col_zoom, ncol = 1),
    ggarrange(plotlist = plt_col_clues, ncol = 1),
    widths = c(2, 3), ncol = 2,
    labels = row_label, font.label = list(color = "grey45"), hjust = 0, vjust = 0
  )
})

# arrange the columns into a grid
plt_grid <- ggarrange(plotlist = plt_columns, ncol = length(focal_regions))

# finally, top the grid with the composite Manhattan plot
main_fig <- ggarrange(
  plotlist = list(plt_man, plt_grid),
  heights = c(1, 4),
  ncol = 1,
  nrow = 2
)

# main_fig

# page dims 183(w) x 247(h) mm
ggsave(filename = argv$output, plot = main_fig, width = 183 * (length(focal_regions) / 3), heigh = 247 * .8, units = "mm", scale = 2, limitsize = FALSE)
