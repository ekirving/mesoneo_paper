#!/usr/bin/env Rscript

# Author:    Evan K. Irving-Pease
# Copyright: Copyright 2020, University of Copenhagen
# Email:     evan.irvingpease@gmail.com
# License:   MIT

quiet <- function(x) {
  suppressMessages(suppressWarnings(x))
}
quiet(library(readr))
quiet(library(dplyr))
quiet(library(stringr))
quiet(library(tidyr))
quiet(library(ggplot2))
quiet(library(bedr))
quiet(library(shadowtext))
quiet(library(ggrepel))
quiet(library(RcppCNPy))
quiet(library(tibble))
quiet(library(zoo))
quiet(library(directlabels))
quiet(library(RcppRoll))

load_clues <- function(clues_tsv, pairs_tsv, unmapped_tsv, flipped_tsv) {

  # all the CLUES results from the GWAS/Control ascertainment
  clues <- read_tsv(clues_tsv, col_types = cols())

  # the GWAS/Control pairings
  pairs <- read_tsv(pairs_tsv, col_types = cols())

  # unmappable and flipped SNPs from Relate
  unmapped <- read_tsv(unmapped_tsv, col_types = cols(pos = "c"))
  flipped <- read_tsv(flipped_tsv, col_types = cols(pos = "c"))

  # calculate p-values from the log-likelihood ratio
  clues$p.value <- pchisq(2 * clues$logLR, df = 1, lower.tail = FALSE)

  # determine if SNPs are GWAS or controls
  clues$type <- ifelse(clues$rsid %in% pairs$gwas, "gwas", "control")

  # set facet ordering
  clues$mode <- factor(clues$mode, levels = c("modern", "ancient"), labels = c("Modern", "Ancient"))
  clues$type <- factor(clues$type, levels = c("gwas", "control"), labels = c("GWAS", "Control"))
  clues$ancestry <- factor(clues$ancestry, levels = c("ALL", "WHG", "EHG", "CHG", "ANA"))

  # drop any modern SNPs that were flipped or unmapped
  clues <- clues %>%
    mutate(pos = paste0(chrom, ":", start)) %>%
    filter(!(mode == "Modern" & (pos %in% unmapped$pos | pos %in% flipped$pos)))

  # filter the data to only retain SNPs for which we have a completed GWAS/Control pair
  clues <- bind_rows(
    rename(pairs, rsid = gwas, pair = neutral),
    rename(pairs, rsid = neutral, pair = gwas)
  ) %>%
    inner_join(clues, by = "rsid") %>%
    select("pair", "mode", "ancestry") %>%
    inner_join(clues, ., by = c("rsid" = "pair", "mode", "ancestry"))

  # tidy up gene name duplicates
  clues <- clues %>%
    rowwise() %>%
    mutate(genes = paste(unique(unlist(str_split(genes, "[,;] *"))), collapse = ","))

  # amazingly, lactase persistence is not in the GWAS catalog!
  clues <- clues %>%
    mutate(gwascat = ifelse(rsid == "rs4988235", paste("Lactase persistence;", gwascat), gwascat))

  clues
}

havester_to_bed <- function(peaks, distance = 1e6) {

  # convert havester output to BED intervals
  bed_peaks <- peaks %>%
    separate(col = range, into = c("start", "end"), sep = "-", convert = TRUE) %>%
    mutate(bed = paste0("chr", chrom, ":", start - 1, "-", end)) %>%
    pull(bed)

  # sort before merging
  bed_peaks <- bedr.sort.region(bed_peaks, verbose = FALSE)

  # merge any regions within +/- 1 Mb of each other
  bed_peaks <- bedr.merge.region(bed_peaks, distance = distance, verbose = FALSE)

  bed_peaks
}

annotate_peaks <- function(data, clues_peaks, mathieson_peaks) {

  # add a BED format column to our df (N.B. BED sorts alphanumerically, so 10 comes before 2)
  data <- data %>%
    arrange(as.character(chrom), start) %>%
    mutate(bed = paste0("chr", chrom, ":", start - 1, "-", end))

  # convert to BED intervals
  mathieson_bed <- havester_to_bed(mathieson_peaks)

  if (nrow(clues_peaks) > 0) {
    clues_bed <- havester_to_bed(clues_peaks)

    # intersect the intervals
    novel_bed <- clues_bed[!clues_bed %in.region% mathieson_bed]
    known_bed <- clues_bed[clues_bed %in.region% mathieson_bed]
    absent_bed <- mathieson_bed[!mathieson_bed %in.region% clues_bed]

    matches <- list()

    if (length(novel_bed)) {
      matches[[1]] <- bedr.join.region(data$bed, novel_bed, verbose = FALSE) %>% mutate(peak = "novel", type = "gwas")
    }

    if (length(known_bed)) {
      matches[[2]] <- bedr.join.region(data$bed, known_bed, verbose = FALSE) %>% mutate(peak = "known", type = "gwas")
    }

    if (length(absent_bed)) {
      matches[[3]] <- bedr.join.region(data$bed, absent_bed, verbose = FALSE) %>% mutate(peak = "absent", type = "gwas")
    }

    # left join the labeled peak regions to the results df
    all_bed <- bind_rows(matches)
  } else {
    clues_bed <- vector()
    novel_bed <- vector()
    known_bed <- vector()
    absent_bed <- mathieson_bed

    all_bed <- bedr.join.region(data$bed, absent_bed, verbose = FALSE) %>% mutate(peak = "absent", type = "gwas")
  }

  all_bed$type <- factor(all_bed$type, levels = c("gwas", "control"), labels = c("GWAS", "Control"))

  data <- all_bed %>%
    filter(V4 != ".") %>%
    mutate(region = paste0(V4, ":", V5, "-", V6)) %>%
    select(index, region, peak, type) %>%
    left_join(data, ., by = c("bed" = "index", "type" = "type"))

  # return the df and the regions
  list("data" = data, "all" = clues_bed, "novel" = novel_bed, "known" = known_bed, "absent" = absent_bed)
}


# helper function to return the nth lowest value
nth_min <- function(x, n) {
  sort(unique(x))[n]
}

label_top_snps <- function(data, facets, num_label, p.threshold) {

  # label the SNPs
  data <- data %>%

    # label the top N SNPs in each peaks and for each facet
    drop_na(region) %>%
    group_by_at(c("region", facets)) %>%
    summarise(cutoff = nth_min(p.value, num_label)) %>%

    # add this info to the initial dataset
    left_join(data, ., by = c("region", facets)) %>%

    # add the annotation flag
    mutate(show_label = (p.value < p.threshold & p.value <= cutoff)) %>%
    select(-cutoff)

  data
}

manhattan_plot <- function(data, facets, p.genomewide, p.bonferroni, num_label = 1, show_lables = TRUE, top_snps = c(), label_colors = list(), highlight_peaks = TRUE, show_count = FALSE, show_legend = FALSE, show_strip = TRUE, wrap = FALSE, composite = FALSE, colors = NA, gene_names = NA, region_name = NA, max_height = NA, plot_title = NA, size_snps = TRUE) {

  # turn off the new summarise warnings
  options(dplyr.summarise.inform = FALSE)

  if (show_lables) {
    if (length(top_snps)) {
      # label just the list of supplied SNPs
      data <- data %>%
        mutate(show_label = (rsid %in% top_snps))
    } else {
      # label the top SNPs in each chromosome (after filtering)
      data <- label_top_snps(data, facets, num_label, p.bonferroni)
    }
  } else {
    data$show_label <- FALSE
  }

  if ( !("merged_peaks" %in% colnames(data)) ) {
      data <- mutate(data, merged_peaks=ifelse(peak != "absent", region, NA))
  }

  data <- data %>%
    # add the labels
    mutate(snp_label = ifelse(show_label, rsid, NA), label_color = replace_na(label_colors[rsid], "black")) %>%

    # size the dots by their `s` coefficient
    mutate(snp_size = ifelse(size_snps && p.value < p.bonferroni & !is.na(merged_peaks), abs(s) * 100, 1)) %>%

    # divide positions by 10, to avoid integer overflow issue with cumulative x-axis position
    mutate(x.pos = start / 10)

  # prepare the dataset
  data <- data %>%

    # compute chromosome size
    group_by(chrom) %>%
    summarise(chr_len = max(x.pos)) %>%

    # calculate cumulative position of each chromosome
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%

    # add this info to the initial dataset
    left_join(data, ., by = "chrom") %>%

    # add a cumulative position of each SNP
    arrange(chrom, x.pos) %>%
    mutate(BPcum = x.pos + tot)

  # center the x-axis labels
  axisdf <- data %>%
    group_by(chrom) %>%
    summarize(center = (max(BPcum) + min(BPcum)) / 2)

  if (wrap) {
    # add a dummy max.height value so we can fix the y-axis across facet_wrap rows
    if (!is.na(max_height)) {
      data$max.height <- max_height
    } else {
      data <- data %>%
        group_by(merged_peaks) %>%
        mutate(max.height = max(-log10(p.value))) %>%
        ungroup()
    }
  }

  # make the plot
  plt <- ggplot(data, aes(x = BPcum, y = -log10(p.value)))

  # subset the points to show peaks and outliers differently
  neutral_snps <- filter(data, p.value >= p.bonferroni & is.na(merged_peaks))
  outlier_snps <- filter(data, p.value < p.bonferroni & is.na(merged_peaks))
  significant_snps <- filter(data, !is.na(merged_peaks))

  if (composite) {

    # draw a box around each of the merged peaks
    peaks <- data %>%
      drop_na(merged_peaks) %>%
      group_by(merged_peaks) %>%
      summarise(xmin = min(BPcum) - 1e6, xmax = max(BPcum) + 1e6, ymax = max(-log10(p.value) + 1)) %>%
      arrange(xmin) %>%
      mutate(xpos = (xmin + xmax) / 2, ypos = ymax + 0.5, label = paste0("Peak ", row_number())) %>%
      # manually dodge peak 13/14 labels
      mutate(xpos = ifelse(row_number() == 13, xpos - 1e6, xpos)) %>%
      mutate(xpos = ifelse(row_number() == 14, xpos + 2e6, xpos), ypos = ifelse(row_number() == 14, ypos - 1, ypos)) %>%
      # add the gene names
      left_join(., gene_names, by = "merged_peaks") %>%
      # only display boxes for GWAS SNPs
      mutate(type = factor("GWAS", levels = c("GWAS", "Control")))

    plt <- plt +
      geom_point(data = neutral_snps, color = "grey", alpha = 0.8) +
      geom_point(data = outlier_snps, color = "grey", alpha = 0.4) +
      geom_rect(data = peaks, inherit.aes = FALSE, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = ymax, group = merged_peaks), alpha = 0.2, fill = "grey") + # , color="grey") +
      geom_point(data = significant_snps, aes(color = ancestry)) +
      geom_shadowtext(data = peaks, aes(x = xpos, y = ypos, label = genes), hjust = "left", angle = 45, color = "black", bg.colour = "white") + # , fontface = "italic") +
      scale_color_manual(values = colors, name = "Ancestry") # , drop = FALSE
  } else {
    if (highlight_peaks) {
      plt <- plt +
        geom_point(data = neutral_snps, color = "grey", alpha = 0.8) +
        geom_point(data = outlier_snps, color = "grey", alpha = 0.4) +
        geom_point(data = significant_snps, aes(color = ancestry), size = significant_snps$snp_size)
    } else {
      plt <- plt +
        geom_point(aes(color = as.factor(chrom)), alpha = 0.8)
    }

    if (wrap) {
      # color the SNPs by their ancestry
      plt <- plt + scale_color_manual(values = colors)
    } else {
      # use an alternating colour scheme
      plt <- plt + scale_color_manual(values = rep(c("grey", "skyblue"), 22))
    }
  }

  if (!is.na(region_name)) {
    # set the region name
    axisdf$chrom <- region_name
  }

  plt <- plt +

    # custom axes
    scale_x_continuous(label = axisdf$chrom, breaks = axisdf$center, expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(add = c(0, ifelse(composite, 5, 1)))) + # add a little padding at the top

    xlab("") +
    geom_hline(yintercept = -log10(p.genomewide), linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(p.bonferroni), linetype = "dashed", color = "grey") +

    # set min height of the y-axis
    # expand_limits(y=-log10(p.genomewide)+0.5) +

    # customise the theme
    theme_bw() +
    theme(
      legend.position = ifelse(show_legend, "right", "none"),
      legend.title = element_text(size = rel(1.1)),
      legend.text = element_text(size = rel(1)),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
    )

  if (show_lables) {
    # add label using ggrepel to avoid overlapping
    plt <- plt +
      geom_text_repel(
        aes(label = snp_label),
        color = data$label_color,
        max.overlaps = Inf,
        point.size = ifelse(length(significant_snps$snp_size) > 0, significant_snps$snp_size, 1),
        size = 3.5,
        min.segment.length = 0,
        box.padding = 1,
        bg.color = "white",
        na.rm = TRUE
      )
  }

  if (!is.na(plot_title)) {
    plt <- plt +
      ggtitle(plot_title) +
      theme(plot.title = element_text(hjust = 0.5, color = "grey45"))
  }

  if (!show_strip) {
    plt <- plt +
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
  }

  if (length(facets) == 1) {
    # display as a wrap
    plt <- plt + facet_wrap(formula(paste("~", facets[1])), ncol = 1)
  } else if (wrap) {
    plt <- plt +
      # include dummy max height
      geom_blank(aes(y = max.height)) +

      # display as a wrap
      facet_wrap(formula(paste(facets[1], "~", facets[2])), scales = "free", ncol = wrap, labeller = labeller(.multi_line = FALSE))
  } else {
    # display as a grid
    plt <- plt + facet_grid(formula(paste(facets[1], "~", facets[2])))
  }

  if (show_count) {
    # get the count of points in each facet
    cell_lables <- data %>%
      group_by_at(facets) %>%
      summarise(label = paste0("n=", n()))

    # add the per-item labels
    plt <- plt + geom_text(size = 3, data = cell_lables, mapping = aes(x = 0, y = Inf, label = label), hjust = 0, vjust = 1.5)
  }

  plt
}

clues_load_data <- function(dataset, population, rsid, mode, ancestry) {
  model_path <- paste0("clues/", dataset, "/", population, "/", rsid, "/", dataset, "-", population, "-", rsid, "-", mode, "-", ancestry, "-any")

  # load the model data
  epochs <- npyLoad(paste0(model_path, ".epochs.npy"))
  freqs <- npyLoad(paste0(model_path, ".freqs.npy"))
  logpost <- npyLoad(paste0(model_path, ".post.npy"), dotranspose = F)

  # add column names
  colnames(logpost) <- paste0("V", seq(ncol(logpost)))

  model <- as_tibble(logpost) %>%

    # convert posterior densities from log-likelihoods
    exp() %>%

    # add the RefSeq ID
    add_column(rsid = rsid, .before = 1) %>%

    # add the frequency labels
    add_column(freq = freqs, .before = 2) %>%

    # add the title heights (and a little padding to make sure there are no gaps)
    # tile heights are not equal because there is higher sampling density near 0 and 1
    add_column(height = diff(c(0, freqs)) + 1e-4, .before = 3) %>%

    # pivot the columns into long format
    pivot_longer(-c(rsid, freq, height), names_to = "epoch", values_to = "density") %>%

    # convert the column names into epochs (and switch the direction of time)
    mutate(epoch = -epochs[as.numeric(gsub("V", "", epoch))])

  model
}

clues_trajectory <- function(dataset, population, rsid, mode, ancestry, ancestral = NULL, smooth = 0) {

  # rsid <- "rs79117855"

  # load the model
  model <- clues_load_data(dataset, population, rsid, mode, ancestry)

  # extract the maximum posterior trajectory
  traj <- model %>%
    group_by(rsid, epoch) %>%
    top_n(1, density) %>%
    ungroup() %>%
    arrange(rsid, epoch)

  # apply a little smoothing to the jagged steps in the trajectory
  if (smooth) {
    traj$freq <- rollapply(c(traj$freq, rep_len(NA, smooth - 1)), width = smooth, by = 1, FUN = mean, na.rm = TRUE, align = "left")
  }

  if (!is.null(ancestral)) {
    metadata <- fromJSON(file = paste0("variants/metadata/GRCh37/", rsid, ".json"))

    # re-polarise, if necessary
    if (metadata$ancestral != ancestral) {
      traj$freq <- 1 - traj$freq
    }
  }

  traj
}

clues_plot <- function(dataset, population, rsid_list, mode, ancestry, snp_colors, geom = "point", title = "", gen_time = 28, max_age = 13665, smooth = 10, ancestral = list()) {

  # dataset <- "ancestral_paths_v3"
  # population <- "all"
  # rsid_list <- c("rs12722987", "rs4988235", "rs79117855")
  # mode <- "ancient"
  # ancestry <- "WHG"
  # title <- "Western Hunter-gatherers"

  # load all the trajectories
  traj <- bind_rows(
    lapply(rsid_list, function(rsid) {
      clues_trajectory(dataset, population, rsid, mode, ancestry, smooth, ancestral = ancestral[[rsid]])
    })
  ) %>% arrange(desc(row_number())) # reverse order, so best is on the top

  # constrain the extent of the plotting
  xmin <- round(-max_age / gen_time)
  xmax <- max(traj$epoch)
  xbreaks <- -seq(-xmax, -xmin, round(2000 / gen_time))
  xlabels <- round(xbreaks * gen_time / 1000)

  # the approximate time (in generations) during which these ancestries existed as discrete populations
  ancestry_epochs <- list(
    "ALL" = c(-150, 0),
    "WHG" = c(-500, -250),
    "EHG" = c(-500, -200),
    "CHG" = c(-800, -200),
    "ANA" = c(-800, -250)
  )

  shade <- unlist(ancestry_epochs[ancestry])

  plt <- traj %>%

    # plot the heatmap
    ggplot(aes(x = epoch, y = freq)) +

    # shade the ancestry epoch
    geom_rect(aes(xmin = max(shade[1], xmin), xmax = min(shade[2], xmax), ymin = 0, ymax = Inf), fill = "#F4F4F4", alpha = 0.1)

  if (geom == "point") {
    plt <- plt +

      # plot the maximum posterior trajectory
      geom_point(aes(color = rsid, alpha = density), cex = 1, na.rm = TRUE) +

      # set the SNP colors
      scale_color_manual(values = snp_colors) +

      # print the labels
      geom_dl(aes(label = rsid, color = rsid), method = list(dl.trans(x = x + 0.1), "last.qp", cex = 0.8), na.rm = TRUE)
  } else if (geom == "line") {
    plt <- plt +

      # plot the maximum posterior trajectory
      geom_line(aes(color = rsid), cex = 0.5)
  }

  plt <- plt +

    # set the axis breaks
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), position = "left", expand = expansion(add = c(0.03, 0))) +
    scale_x_continuous(limits = c(xmin, xmax), breaks = xbreaks, labels = xlabels, expand = expansion(add = c(0, 230))) +

    labs(title = title) +
    ylab("DAF") +
    xlab("kyr BP") +

    # basic styling
    theme_minimal() +
    theme(
      legend.position = "none",
      # panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )

  plt
}

clues_plot_heatmap <- function(prefix, gen_time = 28, max_age = 13665) {

  # load the model data
  epochs <- npyLoad(paste0(prefix, ".epochs.npy"))
  freqs <- npyLoad(paste0(prefix, ".freqs.npy"))
  logpost <- npyLoad(paste0(prefix, ".post.npy"), dotranspose = F)

  # constrain the extent of the plotting
  xmin <- min(epochs)
  xmax <- min(max(epochs), round(max_age / gen_time))
  xbreaks <- seq(xmin, xmax + 1, round(1000 / gen_time))
  xlabels <- round(xbreaks * gen_time / 1000)

  df <- as_tibble(logpost) %>%

    # convert posterior densities from log-likelihoods
    exp() %>%

    # add the frequency labels
    add_column(freq = freqs, .before = 1) %>%

    # add the title heights (and a little padding to make sure there are no gaps)
    # tile heights are not equal because there is higher sampling density near 0 and 1
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

    # labs(title = prefix) +
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

  plt
}
