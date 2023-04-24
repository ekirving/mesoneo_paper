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
quiet(library(english))
quiet(library(jsonlite))
quiet(library(knitr))

# load some helper functions
source("scripts/clues_utils.R")

ancestries <- c("ALL", "WHG", "EHG", "CHG", "ANA")

genes <- read_tsv("clues/ancestral_paths_v3-all-filtered-clues_report-genes.bed", 
                  col_names = c("snp_chr", "snp_start", "snp_end", "rsid", "gene_chr", "gene_start", "gene_end", "closest_gene"), 
                  col_types = cols()) %>%
    select(rsid, closest_gene) %>%
    # handle multiple closest genes
    group_by(rsid) %>%
    summarise(closest_gene = paste0(unique(mixedsort(closest_gene)), collapse = " / "), .groups="drop")

clues <- read_tsv("clues/ancestral_paths_v3-all-filtered-clues_report.tsv", col_types = cols())

# join the closest genes
clues <- clues %>% inner_join(genes, by="rsid")

# get the GWAS and Control pairs
pairs <- read_tsv("variants/ancestral_paths_v3-all-pairs.tsv", col_types = cols())

# set the type field
clues <- clues %>%
    mutate(type = ifelse(rsid %in% pairs$gwas, "gwas", "control"))

# add a BED format column to our df (N.B. BED sorts alphanumerically, so 10 comes before 2)
clues <- clues %>%
    arrange(as.character(chrom), start) %>%
    mutate(bed = paste0("chr", chrom, ":", start - 1, "-", end))

# load all the GWAS Catalog entries
gwascat <- read_tsv("gwascat/gwas_catalog_v1.0.2-associations_e109_r2023-04-07.tsv", guess_max=1e6, quote = "") %>%
    # only keep genome-wide significant
    filter(`P-VALUE` <= 5e-8) %>%
    # select the columns we need
    select(rsid=`SNPS`, phenotype=`DISEASE/TRAIT`, pubmedid=`PUBMEDID`) %>% 
    # split rows where multiple SNPs are reported
    separate_rows(rsid) %>%
    # only keep associations with an actual rsID
    filter(grepl("rs\\d+", rsid)) %>%
    # drop any duplicate entries
    unique()

# group the associations for each rsID
assoc <- gwascat %>%
    group_by(rsid, phenotype) %>% 
    summarise(pubmedid=paste0(sort(pubmedid), collapse = "; "), num_assoc=n()) %>%
    mutate(phenotype=paste0(tolower(substring(phenotype, 1,1)), substring(phenotype, 2))) %>%
    mutate(text=paste0(phenotype, " (", pubmedid, ")")) %>%
    group_by(rsid) %>%
    arrange(desc(num_assoc)) %>%
    summarise(text=combine_words(text))

# write output to this file
sink("supplement-results.txt")

# --------------------------------------------------------------------------------------------------------------
# make the results chapter `Selection in 1000G EUR`
# --------------------------------------------------------------------------------------------------------------
modern <- clues %>% filter(mode=="modern")
modern_all_peaks <- read_tsv("clues/ancestral_paths_v3-all-modern-ALL-filtered-sweeps.tsv", col_types = cols())

# count the SNPs
mod_num_gwas <- modern %>% filter(type=="gwas") %>% pull(rsid) %>% unique() %>% length()
mod_num_control <- modern %>% filter(type=="control") %>% pull(rsid) %>% unique() %>% length()

# count the significant SNPs
mod_num_signif <- modern %>% filter(p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()
mod_num_gwas_signif <- modern %>% filter(type=="gwas" & p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()
mod_num_control_signif <- modern %>% filter(type=="control" & p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()

# count the sweeps
mod_num_sweeps <- modern_all_peaks %>% nrow()
mod_num_gwas_sweeps <- modern_all_peaks %>% filter(type=="gwas") %>% nrow()
mod_num_control_sweeps <- modern_all_peaks %>% filter(type=="control") %>% nrow()

mod_paired <- pairs %>%
    inner_join(select(modern, c(rsid, p.value)), by = c("gwas" = "rsid")) %>%
    inner_join(select(modern, c(rsid, p.value)), by = c("neutral" = "rsid"), suffix = c(".gwas", ".neutral"))

mod_wilcox <- wilcox.test(mod_paired$p.value.gwas, mod_paired$p.value.neutral, paired = TRUE, alternative = "less")

# get the most significant SNP
mod_top_snp <- modern %>% filter(type=="gwas") %>% group_by() %>% slice_min(p.value, with_ties = FALSE)

# join all the associations, including any for synonymous rsIDs
mod_top_snp <- mod_top_snp %>%
    left_join(
        bind_rows(
            lapply(mod_top_snp$rsid, function(id) {
                synonyms <- c(id, fromJSON(paste0("variants/ensembl/GRCh37/var/", id, ".json"))$synonyms)
                assoc %>% filter(rsid %in% synonyms) %>% mutate(rsid=id)
            })
        ) %>% group_by(rsid) %>% summarise(text=paste0(text, collapse = ", ")),
        by = "rsid"
    ) %>% mutate(severity=str_replace_all(severity, "_", " "))

# compose the summary paragraph
mod_intro <- paste0(
    "CLUES analysis of all GWAS (n=", sprintf("%d", mod_num_gwas),") and Control group SNPs (n=", mod_num_control, 
    "), which passed all quality control filters, in the 1000G Project populations FIN, GBR, and TSI, ",
    "identified ", mod_num_signif," genome-wide significant SNPs (p<5e-8); ", mod_num_gwas_signif," in the GWAS group (", 
    sprintf("%.2f", (mod_num_gwas_signif/mod_num_signif)*100, 4) ,"%) and ", mod_num_control_signif," in the Control group (", sprintf("%.2f", (mod_num_control_signif/mod_num_signif)*100, 4) ,
    "%). Within the GWAS group, we identified ", mod_num_gwas_sweeps, " genome-wide significant sweep regions, ",
    "and ", mod_num_control_sweeps," in the Control group (see Fig. S2a.2; Supplementary Table S2a.XX). ",
    "Despite the general lack of genome-wide significant SNPs, the GWAS group was significantly enriched for evidence of selection when compared to the Control group ",
    "(Wilcoxon signed-rank test, p-value < ", sprintf("%.2e", mod_wilcox$p.value), ")."
)

lct_snp <- modern %>% filter(rsid=="rs4988235")

mod_summary <- mod_top_snp %>%
    mutate(text=paste0(
        "The most significant SNP was ", rsid, " (", closest_gene, "; p=", sprintf("%.2e", p.value) ,"; s=", sprintf("%.3f", s), "), an ", 
        severity, " associated with ", text," in the GWAS Catalog (r2023-04-07). Other genome-wide significant SNPs within the surrounding region ",
        "include the lactase persistence SNP rs4988235 (p=", sprintf("%.2e", lct_snp$p.value), "; s=", sprintf("%.3f", lct_snp$s), "), ",
        "which has been widely reported as a target of strong selection in West Eurasians (Enattah et al. 2002; Bersaglieri et al. 2004)"
    ))


mod_manhat_caption <- "Figure S2a.2. Manhattan plot of the p-values from running CLUES on an ARG containing all samples in FIN, GBR, and TSI from 1000G Phase 3, for (a) GWAS SNPs from the GWAS Catalog; and (b) Control SNPs, frequency paired with the GWAS SNPs."

cat("Selection in 1000G EUR", "\n\n")
cat(mod_intro, "\n\n")
cat(mod_summary$text, "\n\n")
cat(mod_manhat_caption, "\n\n")

# --------------------------------------------------------------------------------------------------------------
# make the results chapter `Selection in aDNA Time Series`
# --------------------------------------------------------------------------------------------------------------

ancient <- clues %>% filter(mode=="ancient" & ancestry=="ALL")
ancient_all_peaks <- read_tsv("clues/ancestral_paths_v3-all-ancient-ALL-filtered-sweeps.tsv", col_types = cols())
ancient_all_peaks_bed <- sweep_to_bed(ancient_all_peaks)

# left join the sweep regions
ancient <- bedr.join.region(ancient$bed, ancient_all_peaks_bed, verbose = FALSE) %>%
    filter(V4 != ".") %>%
    mutate(region = paste0(V4, ":", V5, "-", V6)) %>%
    select(index, region) %>%
    left_join(ancient, ., by = c("bed" = "index"))

# count the SNPs
anc_num_gwas <- ancient %>% filter(type=="gwas") %>% pull(rsid) %>% unique() %>% length()
anc_num_control <- ancient %>% filter(type=="control") %>% pull(rsid) %>% unique() %>% length()

# count the significant SNPs
anc_num_signif <- ancient %>% filter(p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()
anc_num_gwas_signif <- ancient %>% filter(type=="gwas" & p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()
anc_num_control_signif <- ancient %>% filter(type=="control" & p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()

# count the sweeps
anc_num_sweeps <- ancient_all_peaks %>% nrow()
anc_num_gwas_sweeps <- ancient_all_peaks %>% filter(type=="gwas") %>% nrow()
anc_num_control_sweeps <- ancient_all_peaks %>% filter(type=="control") %>% nrow()

# summarise the sweeps
ancient_summary <- ancient %>%
    drop_na(region) %>% 
    group_by(region) %>% 
    slice_min(p.value, with_ties = FALSE) %>% 
    ungroup() %>% 
    arrange(chrom, start) %>% 
    mutate(peak=row_number()) %>%
    mutate(severity=str_replace_all(severity, "_", " ")) %>%
    select(type, peak, region, closest_gene, rsid, mode, ancestry, p.value, s, severity)

# join all the associations, including any for synonymous rsIDs
ancient_summary <- ancient_summary %>%
    left_join(
        bind_rows(
            lapply(ancient_summary$rsid, function(id) {
                synonyms <- c(id, fromJSON(paste0("variants/ensembl/GRCh37/var/", id, ".json"))$synonyms)
                assoc %>% filter(rsid %in% synonyms) %>% mutate(rsid=id)
            })
        ) %>% group_by(rsid) %>% summarise(text=paste0(text, collapse = ", ")),
        by = "rsid"
    )

ancient_content <- ancient_summary %>%
    mutate(title=paste0("Peak ", peak, ": ", closest_gene)) %>%
    mutate(caption1=paste0("Figure S2a.", 4+(peak-1)*2, ". CLUES plot of the aDNA time series analysis for all West Eurasian samples in the imputed dataset, showing the posterior probability of the derived allele frequency trajectory for ", rsid, " the most significant SNP in the selection peak spanning ", region, ".")) %>%
    mutate(caption2=paste0("Figure S2a.", 5+(peak-1)*2, ". CLUES plot of the aDNA time series analysis for all West Eurasian samples in the genotype likelihood dataset (i.e., non-imputed), showing the posterior probability of the derived allele frequency trajectory for ", rsid, ".")) %>%
    mutate(description=paste0("The ", ordinal(peak), " peak spanned the region ", region, ", with the most significant SNP being ", rsid, " (", closest_gene, "; p=", sprintf("%.2e", p.value) ,"; s=", sprintf("%.3f", s), "), an ", severity, " associated with ", text," in the GWAS Catalog (r2023-04-07).")) %>%
    mutate(content=paste(title, caption1, caption2, description, sep = "\n\n")) %>%
    pull(content)

# compose the summary paragraph
anc_intro <- paste0(
    "CLUES analysis of all GWAS (n=", sprintf("%d", anc_num_gwas),") and Control group SNPs (n=", anc_num_control, 
    "), which passed all quality control filters in the aDNA time-series dataset, ",
    "identified ", anc_num_signif," genome-wide significant SNPs (p<5e-8); ", anc_num_gwas_signif," in the GWAS group (", 
    sprintf("%.2f", (anc_num_gwas_signif/anc_num_signif)*100, 4) ,"%) and ", anc_num_control_signif," in the Control group (", sprintf("%.2f", (anc_num_control_signif/anc_num_signif)*100, 4) ,
    "%). Within the GWAS group, we identified ", anc_num_gwas_sweeps, " genome-wide significant sweep regions, ",
    "and ", anc_num_control_sweeps," in the Control group (see Fig. S2a.3; Supplementary Table S2a.XX). "
)

anc_manhat_caption <- "Figure S2a.3. Manhattan plot of the p-values from running CLUES on an aDNA time series from all West Eurasian samples in the imputed dataset, for (a) GWAS SNPs from the GWAS Catalog; and (b) Control SNPs, frequency paired with the GWAS SNPs."

cat("Selection in aDNA Time Series", "\n\n")
cat(anc_intro, "\n\n")
cat(anc_manhat_caption, "\n\n")
cat(paste0(ancient_content, collalse="\n\n"))

# --------------------------------------------------------------------------------------------------------------
# make the results chapter `Selection in simulations with Ancestral Paintings`
# --------------------------------------------------------------------------------------------------------------

simulated <- read_tsv("clues/simulated_relate_painted-all-simulated-clues_report.tsv", col_types = cols())
simulated_peaks <- read_tsv("clues/simulated_relate_painted-all-simulated-filtered-sweeps-merged.tsv", col_types = cols())

# count the SNPs
sim_num_snps <- simulated %>% pull(rsid) %>% unique() %>% length()

# count the significant SNPs
sim_num_signif <- simulated %>% filter(p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()

# count the sweeps
sim_num_sweeps <- simulated_peaks %>% nrow()

# compose the summary paragraph
sim_intro <- paste0(
    "CLUES analysis of all simulated SNPs (n=", sim_num_snps, "), in the frequency-paired aDNA simulation, identified ", 
    sim_num_signif," (", sprintf("%.2f", (sim_num_signif/sim_num_snps)*100, 4) ,"%) genome-wide significant SNPs (p<5e-8). ",
    "We identified ", sim_num_sweeps, " genome-wide significant sweep regions across all ancestries (see Fig. S2a.26; Supplementary Table S2a.XX), ",
    "indicating that the false positive rate of our ancestry stratified selection analysis is low."
)


sim_manhat_caption <- "Figure S2a.26. Manhattan plot of the p-values from running CLUES on a neutral simulation of chr3, using the inferred ancestral paths of each ancestry. The first row shows results for all ancient samples considered in aggregate, and each subsequent row shows the results conditional on one of the four specific ancestral paintings:  WHG (Western Hunter-gatherers), EHG (Eastern Hunter-gatherers), CHG (Caucasus Hunter-gatherers), and ANA (Anatolian Farmers)."

cat("Selection in simulations with Ancestral Paintings", "\n\n")
cat(sim_intro, "\n\n")
cat(sim_manhat_caption, "\n\n")


# --------------------------------------------------------------------------------------------------------------
# make the results chapter `Selection in aDNA with Ancestral Paintings`
# --------------------------------------------------------------------------------------------------------------

ancestral <- clues %>% filter(mode=="ancient")
merged_peaks <- read_tsv("clues/ancestral_paths_v3-all-ancient-filtered-sweeps-merged.tsv", col_types = cols())
merged_peaks_bed <- sweep_to_bed(merged_peaks)

# left join the sweep regions
ancestral <- bedr.join.region(ancestral$bed, merged_peaks_bed, verbose = FALSE) %>%
    filter(V4 != ".") %>%
    mutate(region = paste0(V4, ":", V5, "-", V6)) %>%
    select(index, region) %>%
    left_join(ancestral, ., by = c("bed" = "index")) %>%
    unique()

# count the SNPs
num_gwas <- ancestral %>% filter(type=="gwas") %>% pull(rsid) %>% unique() %>% length()
num_control <- ancestral %>% filter(type=="control") %>% pull(rsid) %>% unique() %>% length()

# count the significant SNPs
num_signif <- ancestral %>% filter(p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()
num_gwas_signif <- ancestral %>% filter(type=="gwas" & p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()
num_control_signif <- ancestral %>% filter(type=="control" & p.value <= 5e-8) %>% pull(rsid) %>% unique() %>% length()

# count the sweeps
num_sweeps <- merged_peaks %>% nrow()
num_gwas_sweeps <- merged_peaks %>% filter(type=="gwas") %>% nrow()
num_control_sweeps <- merged_peaks %>% filter(type=="control") %>% nrow()

NUM_TOP_SNPS <- 5

# get the top SNPs across all the marginal ancestries
top_marginal_snps <- ancestral %>%
    # drop any SNPs that are not in the actual peak
    drop_na(region) %>%
    # only keep significant SNPs
    filter(p.value <= 5e-8) %>%
    # get the top SNPs in each marginal ancestry
    group_by(region, ancestry) %>%
    slice_min(p.value, n = NUM_TOP_SNPS, with_ties = FALSE) %>%
    mutate(ancestry_rank = row_number()) %>%
    ungroup() %>%
    # order by the ancestry rank and the marginal p-value
    arrange(region, ancestry_rank, p.value) %>%
    # drop any duplicates across ancestries
    select(region, rsid) %>%
    unique() %>%
    # get the top SNPs across all ancestries
    group_by(region) %>%
    mutate(peak_rank = row_number()) %>%
    filter(peak_rank <= NUM_TOP_SNPS) %>%
    # get the most significant ancestry for each of the top SNPs
    inner_join(ancestral, by = c("region", "rsid")) %>%
    group_by(rsid) %>%
    slice_min(p.value, with_ties = FALSE) %>%
    select(region, rsid, closest_gene, severity, peak_rank, ancestry, p.value, s) %>%
    # order by the peak rank
    arrange(region, peak_rank)

# join all the associations, including any for synonymous rsIDs
top_marginal_snps <- top_marginal_snps %>%
    left_join(
        bind_rows(
            lapply(top_marginal_snps$rsid, function(id) {
                synonyms <- c(id, fromJSON(paste0("variants/ensembl/GRCh37/var/", id, ".json"))$synonyms)
                assoc %>% filter(rsid %in% synonyms) %>% mutate(rsid=id)
            })
        ) %>% group_by(rsid) %>% summarise(text=paste0(text, collapse = ", ")),
        by = "rsid"
    )

# summarise each SNP
snps_summary <- top_marginal_snps %>%
    mutate(severity=str_replace_all(severity, "_", " ")) %>%
    mutate(text=paste0(rsid, " (", closest_gene, "; ", ancestry, ": p=", sprintf("%.2e", p.value) ,"; s=", sprintf("%.3f", s), "), an ", severity, " associated with ", text," in the GWAS Catalog (r2023-04-07).")) %>%
    group_by(region) %>%
    summarise(text=paste(text, collapse = "\n"))

# summarise the sweeps
ancestral_summary <- ancestral %>%
    drop_na(region) %>% 
    group_by(region) %>% 
    slice_min(p.value, with_ties = FALSE) %>% 
    ungroup() %>% 
    arrange(chrom, start) %>% 
    mutate(peak=row_number()) %>%
    select(type, peak, region, closest_gene) %>%
    left_join(snps_summary, by="region")

ancestal_content <- ancestral_summary %>%
    mutate(title=paste0("Peak ", peak, ": ", closest_gene)) %>%
    mutate(caption=paste0("Figure S2a.", 29+peak, ". Selection at the ", closest_gene, " locus, spanning ", region, ". Results for the pan-ancestry analysis (ALL) plus each of the four marginal ancestries: Western hunter-gatherers (WHG), Eastern hunter-gatherers (EHG), Caucasus hunter-gatherers (CHG) and Anatolian farmers (ANA). Row one shows zoomed Manhattan plots of the p-values for each ancestry, and row two shows allele trajectories for the top SNPs across all ancestries.")) %>%
    mutate(description=paste0("The ", ordinal(peak), " peak spanned the region ", region, ", with the five most significant SNPs across ancestries: ", "\n", text)) %>%
    mutate(content=paste(title, caption, description, sep = "\n\n")) %>%
    pull(content)

# compose the summary paragraph
ancestral_intro <- paste0(
    "CLUES analysis of all GWAS (n=", sprintf("%d", num_gwas),") and Control group SNPs (n=", num_control, ") in the aDNA with Ancestral Paintings ",
    "dataset identified ", num_signif," genome-wide significant SNPs (p<5e-8); ", num_gwas_signif," in the GWAS group (", 
    sprintf("%.2f", (num_gwas_signif/num_signif)*100, 4) ,"%) and ", num_control_signif," in the Control group (", sprintf("%.2f", (num_control_signif/num_signif)*100, 4) ,
    "%). Within the GWAS group, we identified ", num_gwas_sweeps, " non-overlapping genome-wide significant sweep regions across all ancestries, ",
    "and ", num_control_sweeps," in the Control group (see Fig. S2a.27 and S2a.28; Supplementary Table S2a.XX)."
)

ancestral_gwas_caption <- "Figure S2a.27. Manhattan plot of the p-values from running CLUES on an aDNA time series conditioned on ancestry paintings from all West Eurasian samples in the imputed dataset for GWAS SNPs from the GWAS Catalog. The first row shows results for all ancient samples considered in aggregate, and each subsequent row shows the results conditional on one of the four specific ancestral paintings: WHG (Western Hunter-gatherers), EHG (Eastern Hunter-gatherers), CHG (Caucasus Hunter-gatherers), and ANA (Anatolian Farmers)."
ancestral_control_caption <- "Figure S2a.28. Manhattan plot of the p-values from running CLUES on an aDNA time series conditioned on ancestry paintings from all West Eurasian samples in the imputed dataset for Control SNPs, frequency paired with the GWAS SNPs. The first row shows results for all ancient samples considered in aggregate, and each subsequent row shows the results conditional on one of the four specific ancestral paintings: WHG (Western Hunter-gatherers), EHG (Eastern Hunter-gatherers), CHG (Caucasus Hunter-gatherers), and ANA (Anatolian Farmers)."

cat("Selection in aDNA with Ancestral Paintings", "\n\n")
cat(ancestral_intro, "\n\n")
cat(ancestral_gwas_caption, "\n\n")
cat(ancestral_control_caption, "\n\n")
cat(paste0(ancestal_content, collalse="\n\n"))

# close the sink
sink()

# --------------------------------------------------------------------------------------------------------------
# generate come copy commands to organize the figures for the supplement
# --------------------------------------------------------------------------------------------------------------

cp_ancient <- ancient_summary %>%
    mutate(impute=paste0("cp clues/ancestral_paths_v3/all/", rsid, "/ancestral_paths_v3-all-", rsid, "-ancient-ALL-any.png figs/supplement/ancient/peak_", peak, "_ancestral_paths_v3-all-", rsid, "-ancient-ALL-any.png")) %>%
    mutate(likeli=paste0("cp clues/neo_likelihoods/all/", rsid, "/neo_likelihoods-all-", rsid, "-ancient-ALL-any.png figs/supplement/ancient/peak_", peak, "_neo_likelihoods-all-", rsid, "-ancient-ALL-any.png"))

cat(paste0(cp_ancient %>% pull(impute), collalse="\n"))
cat("\n")
cat(paste0(cp_ancient %>% pull(likeli), collalse="\n"))

# generate the supplementary tables

# Supplementary Table S2a.1 - Modern CLUES
modern %>% select(-c(genes, bed)) %>% write_tsv("Supplementary_Table_S2a.1-Modern_CLUES.tsv")

# Supplementary Table S2a.2 - Simulated CLUES
simulated %>% write_tsv("Supplementary_Table_S2a.2-Simulated_CLUES.tsv")

# Supplementary Table S2a.3 - Pan-ancestry CLUES
ancient %>% select(-c(genes, bed)) %>% rename(sweep=region) %>% write_tsv("Supplementary_Table_S2a.3-Pan-ancestry_CLUES.tsv")

# Supplementary Table S2a.4 - Ancestry stratified CLUES
ancestral %>% select(-c(genes, bed)) %>% rename(sweep=region) %>% write_tsv("Supplementary_Table_S2a.4-Ancestry_stratified_CLUES.tsv")

# Supplementary Table S2a.5 - LCT sweep


