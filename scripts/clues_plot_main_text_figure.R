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
p <- arg_parser("Plot the CLUES main text figure")
p <- add_argument(p, "--data", help = "CLUES report", default = "clues/ancestral_paths_v3-all-filtered-clues_report.tsv")
p <- add_argument(p, "--genes", help = "Nearest gene found in Ensembl", default = "clues/ancestral_paths_v3-all-filtered-clues_report-genes.bed")
p <- add_argument(p, "--pairs", help = "GWAS and neutral SNP pairings", default = "variants/ancestral_paths_v3-all-pairs.tsv")
p <- add_argument(p, "--unmapped", help = "List of modern SNPs unmappable by Relate", default = "relate/1000G_phase3-FIN_GBR_TSI-popsize-allsnps_unmapped.tsv.gz")
p <- add_argument(p, "--flipped", help = "List of modern SNPs flipped by Relate", default = "relate/1000G_phase3-FIN_GBR_TSI-popsize-allsnps_flipped.tsv.gz")
p <- add_argument(p, "--peaks-all", help = "CLUES selection peaks from the pan-ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-ALL-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-whg", help = "CLUES selection peaks from the WHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-WHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-ehg", help = "CLUES selection peaks from the EHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-EHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-chg", help = "CLUES selection peaks from the CHG ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-CHG-filtered-sweeps.tsv")
p <- add_argument(p, "--peaks-ana", help = "CLUES selection peaks from the ANA ancestry analysis", default = "clues/ancestral_paths_v3-all-ancient-ANA-filtered-sweeps.tsv")
p <- add_argument(p, "--merged", help = "Merged selection peaks across ancestries", default = "clues/ancestral_paths_v3-all-ancient-filtered-sweeps-merged.tsv")
p <- add_argument(p, "--mathieson", help = "Selection peaks from Mathieson et al. 2015", default = "mathieson/mathieson-sweeps.tsv")
p <- add_argument(p, "--output", help = "Output file", default = "figs/ancestral_paths_v3-all-filtered-main_figure.png")

argv <- parse_args(p)

# load some helper functions
source("scripts/clues_utils.R")

ancestries <- c("ALL", "WHG", "EHG", "CHG", "ANA")

# the width of the flanking region shown in the small Manhattan plots
REGION_PADDING <- 1e6

# load the nearest gene for each SNP (because the control SNPs lack gene annotations)
ensembl_genes <- read_tsv(argv$genes, col_names = c("snp_chr", "snp_start", "snp_end", "rsid", "gene_chr", "gene_start", "gene_end", "closest_gene"), col_types = cols()) %>%
    select(rsid, closest_gene) %>%
    # handle multiple closest genes
    group_by(rsid) %>%
    summarise(closest_gene = paste0(unique(mixedsort(closest_gene)), collapse = " / "), .groups="drop")

# load all the CLUES results
clues <- load_clues(argv$data, argv$pairs, argv$unmapped, argv$flipped) %>%
    left_join(ensembl_genes, by="rsid")

# load all the CLUES peaks
ancient_all_peaks <- read_tsv(argv$peaks_all, col_types = cols())
ancient_whg_peaks <- read_tsv(argv$peaks_whg, col_types = cols())
ancient_ehg_peaks <- read_tsv(argv$peaks_ehg, col_types = cols())
ancient_chg_peaks <- read_tsv(argv$peaks_chg, col_types = cols())
ancient_ana_peaks <- read_tsv(argv$peaks_ana, col_types = cols())
ancient_merged_peaks <- read_tsv(argv$merged, col_types = cols())

# peaks for the reported p-values in Mathieson_et_al_2015
mathieson_peaks <- read_tsv(argv$mathieson, col_types = cols())

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

all_regions <- tibble(
    merged_peaks = sweep_to_bed(ancient_merged_peaks),
    # add some buffer either side (for display purposes only; does not effect the reported intervals)
    buffered_peaks = sweep_to_bed(ancient_merged_peaks %>% mutate(
            start = ifelse(start > REGION_PADDING, start - REGION_PADDING, 0),
            end = end + REGION_PADDING
        ))
    )

left_join_regions <- function(peaks_joined, regions, column) {
    # join the merged regions
    bedr.join.region(peaks_joined$bed, regions, verbose = FALSE) %>%
        filter(V4 != ".") %>%
        mutate(!!column := paste0(V4, ":", V5, "-", V6)) %>%
        select(index, !!column) %>%
        unique() %>%
        left_join(peaks_joined, ., by = c("bed" = "index"))
}

# join the merged peaks and the buffered peaks
peaks_joined <- left_join_regions(peaks_joined, all_regions$merged_peaks, "merged_peaks")
peaks_joined <- left_join_regions(peaks_joined, all_regions$buffered_peaks, "buffered_peaks")

# for SNPs inside selection peaks, only retain the most significant marginal ancestry
ancestral_min <- peaks_joined %>%
  filter((ancestry == "ALL" & is.na(merged_peaks)) | (ancestry != "ALL" & !is.na(merged_peaks))) %>%
  group_by(rsid) %>%
  slice_min(p.value, n = 1, with_ties = FALSE)

# get the gene names for each peak
gene_names <- peaks_joined %>%
    drop_na(merged_peaks) %>%
    group_by(merged_peaks) %>%
    slice_min(p.value, with_ties = FALSE) %>%
    arrange(chrom, start) %>%
    select(chrom, merged_peaks, buffered_peaks, genes=closest_gene, p.value) %>%
    mutate(genes = ifelse(grepl("chr6:.+", merged_peaks), "HLA", genes))

# use a genome-wide significance threshold
p.genomewide <- 5e-8

# plot the ancestry stratified selection peaks
plt_man <- manhattan_plot(
    ancestral_min, 
    c("mode"), 
    p.genomewide, 
    show_lables = FALSE, 
    show_legend = TRUE, 
    show_strip = FALSE, 
    composite = TRUE, 
    colors = ancestry_colors, 
    gene_names = gene_names,
    size_snps = FALSE
  ) +
  theme(plot.title = element_text(hjust = 0.5))

# plt_man

# ------------------------------------------------------------------------------------------------------------------------------------------------------

# pick three regions for the zoomed-in columns
focal_regions <- gene_names %>%
  filter(genes %in% c("MCM6", "SLC45A2", "FADS2")) %>%
  pull(buffered_peaks) %>%
  mixedsort()

mcm6_region <- gene_names %>% filter(genes == "MCM6") %>% pull(buffered_peaks)

# get the top SNP from each of the marginal ancestries
top_marginal_snps <- peaks_joined %>%
  drop_na(merged_peaks) %>%
  filter(p.value < p.genomewide) %>%
  group_by(buffered_peaks, ancestry) %>%
  slice_min(p.value, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(buffered_peaks, rsid) %>%
  unique()

plt_columns <- lapply(focal_regions, function(focal_region) {
  # focal_region <- focal_regions[[3]]

  # get all the SNPs in this region
  region_data <- peaks_joined %>% filter(buffered_peaks == focal_region)

  # get the list of top SNPs in the region
  region_snps <- top_marginal_snps %>% filter(buffered_peaks == focal_region) %>% pull(rsid)

  if (focal_region == mcm6_region) {
    # force display of the microRNA SNP
    if (!("rs1438307" %in% region_snps)) {
       region_snps <- c(region_snps, "rs1438307")
    }
  }

  # pair the colours with the ordered top SNPs
  label_colors <- setNames(snp_colors[1:length(region_snps)], region_snps)

  # make region name human readable
  region_name <- unlist(strsplit(focal_region, split = "[:-]"))
  region_name <- sprintf("%s:%.1f-%.1f Mb", region_name[1], as.integer(region_name[2]) / 1e6, as.integer(region_name[3]) / 1e6)

  # determine the max height, to keep the y-axis in sync
  max_height <- max(-log10(region_data$p.value))

  # make the zoomed in Manhattan plots
  plt_col_zoom <- lapply(ancestries, function(focal_ancestry) {
    # focal_ancestry <- ancestries[[1]]

    manhattan_plot(
        filter(region_data, ancestry == focal_ancestry),
        c("buffered_peaks", "ancestry"),
        p.genomewide,
        top_snps = region_snps,
        show_strip = FALSE,
        wrap = 1,
        composite = FALSE,
        colors = ancestry_colors,
        label_colors = label_colors,
        region_name = region_name,
        max_height = max_height
    )
  })

  # make the per-ancestry CLUES plots
  plt_col_clues <- lapply(ancestries, function(focal_ancestry) {
    # focal_ancestry <- ancestries[[1]]

    # sort the SNPs by their marginal p-value
    ordered_snps <- region_data %>%
      filter(rsid %in% region_snps & ancestry == focal_ancestry) %>%
      arrange(p.value) %>%
      pull(rsid)

    clues_plot("ancestral_paths_v3", "all", ordered_snps, "ancient", focal_ancestry, label_colors)
  })

  row_label <- filter(gene_names, buffered_peaks == focal_region)$genes

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
