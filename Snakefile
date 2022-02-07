#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import math
import re

import pandas as pd
import pysam
from snakemake.io import expand


configfile: "config.yaml"


include: "rules/tgp.smk"
include: "rules/gwascat.smk"
include: "rules/ensembl.smk"
include: "rules/variants.smk"
include: "rules/bedtools.smk"
include: "rules/reports.smk"
include: "rules/neutral.smk"
include: "rules/simulated.smk"
include: "rules/inv17_h1h2.smk"
include: "rules/relate.smk"
include: "rules/clues.smk"
include: "rules/mathieson.smk"


CLUES_MODES = ["ancient", "modern", "both"]
ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]
AGE_MODEL = ["model", "noage"]
PATHS = ["true", "inferred"]


wildcard_constraints:
    reference="[\w.]+",
    dataset="[\w.]+",
    panel="[\w.]+",
    population="[a-zA-Z0-9.]+",
    paths="|".join(PATHS),
    pops="[A-Z_]+",
    rsid="(rs\d+|chr\d+:\d+)",
    age="|".join(AGE_MODEL),
    mode="|".join(CLUES_MODES),
    ancestry="|".join(ANCESTRIES),
    n="\d+",
    chr="\d+",
    bias="[\w_]+",


def get_paired_snps(dataset, population, snp_type, chrom=None):
    simulated = ["chr3_true_paths", "chr3_inferred_paths"]

    if dataset in simulated:
        # noinspection PyUnresolvedReferences
        pairs = checkpoints.simulated_pair_snps.get()
    else:
        # noinspection PyUnresolvedReferences
        pairs = checkpoints.neutral_pair_snps.get(dataset=dataset, population=population)

    data = pd.read_table(pairs.output[0])

    if chrom is not None:
        # apply the chromosome filter
        data = data[data["chr"] == chrom]

    return data[snp_type].tolist()


def run_all_ancestries(_):
    dataset = config.get("dataset", "ancestral_paths_new")
    population = config.get("population", "all")
    ancestry = config.get("ancestry", "ALL")
    snp_type = config.get("type", "gwas")

    # get all the SNPs
    rsid = get_paired_snps(dataset, population, snp_type, config.get("chr", None))

    return expand(
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-{ancestry}-any.png",
        dataset=dataset,
        population=population,
        rsid=rsid,
        ancestry=ancestry,
    )


def get_unmapped():
    """Get the list of SNPs that Relate could not map"""

    # noinspection PyUnresolvedReferences
    unmapped_rule = checkpoints.relate_is_not_mapping.get(panel="1000G_phase3", pops="FIN_GBR_TSI")
    unmapped = pd.read_table(unmapped_rule.output[0])

    return set(unmapped["rsid"])


def run_all_moderns(_):
    dataset = config.get("dataset", "ancestral_paths_new")
    population = config.get("population", "all")
    snp_type = config.get("type", "gwas")

    # get all the SNPs
    rsid = get_paired_snps(dataset, population, snp_type, config.get("chr", None))

    # drop all the non-mapping SNPs
    rsid = set(rsid) - get_unmapped()

    return expand(
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-modern-ALL-any.png",
        dataset=dataset,
        population=population,
        rsid=rsid,
    )


def run_all_andres(_):
    dataset = "imputed_unfiltered"
    population = "all"

    # load the list of SNPs supplied by Andrés
    data = pd.read_table(config["andres"]["targets"])
    data = data[data["chrom"] != "chrX"]

    rsid = set(data["rsid"].tolist())

    return expand(
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-ALL-any.png",
        dataset=dataset,
        population=population,
        rsid=rsid,
    )


def run_all_inv17_h1h2(_):
    dataset = config.get("dataset", "ancestral_paths_new")
    population = config.get("population", "all")
    ancestries = ANCESTRIES if dataset == "ancestral_paths_new" else "ALL"

    # load the list of SNPs supplied by Andrés
    data = pd.read_table(config["inv17_h1h2"]["targets"])

    # noinspection PyUnresolvedReferences
    info = checkpoints.bcftools_info_info.get(dataset=dataset, population=population)
    info = pd.read_table(info.output[0])

    # drop the non-rsID part of the ID
    info["rsid"] = info["rsid"].str.extract(pat=r"(rs\d+)")

    # filter out the missing SNPs
    data = data[data.SNP.isin(info["rsid"])]
    rsid = set(data["SNP"].tolist())

    return expand(
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-{ancestry}-any.png",
        dataset=dataset,
        population=population,
        rsid=rsid,
        ancestry=ancestries,
    )


def run_all_mathieson(_):
    dataset = "ancestral_paths_new"
    population = "all"

    # noinspection PyUnresolvedReferences
    mathieson = checkpoints.mathieson_significant.get()

    data = pd.read_table(mathieson.output[0])
    all_snps = set(data["ID"].tolist())

    # drop all the non-mapping SNPs
    modern_snps = all_snps - get_unmapped()

    modern = expand(
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-modern-ALL-any.png",
        dataset=dataset,
        population=population,
        rsid=modern_snps,
    )

    ancient = expand(
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-{ancestry}-any.png",
        dataset=dataset,
        population=population,
        rsid=all_snps,
        ancestry=ANCESTRIES,
    )

    return modern + ancient


def run_all_locus(_):
    dataset = config.get("dataset", "ancestral_paths_new")
    population = config.get("population", "all")
    locus = config.get("locus", "chr17:25869029-48635050")
    batch = config.get("batch", "1")

    vcf_file = config["samples"][dataset]["genotypes"]

    # fetch all the SNPs in this locus
    vcf = pysam.VariantFile(vcf_file)

    rsid = []

    chrom, start, end = re.split("[:-]", locus)

    # fetch the rsIDs of all the SNPs in this locus
    for rec in vcf.fetch(chrom.replace("chr", ""), int(start) - 1, int(end)):
        match = re.search("(rs\d+)", rec.id)
        if match:
            rsid.append(match.group(1))

    # apply the batching
    chunk = math.ceil(len(rsid) / 100)
    start = (batch - 1) * chunk

    return expand(
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-{ancestry}-any.png",
        dataset=dataset,
        population=population,
        rsid=rsid[start : start + chunk],
        ancestry=ANCESTRIES,
    )


rule all:
    # produce all the selection reports and pick peaks with manhattan harvester
    input:
        "figs/ancestral_paths_new-all-filtered-main_figure.png",


rule ancestries:
    input:
        run_all_ancestries,


rule modern:
    input:
        run_all_moderns,


rule andres:
    input:
        run_all_andres,


rule inv17_h1h2:
    input:
        run_all_inv17_h1h2,


rule mathieson:
    input:
        run_all_mathieson,


rule locus:
    input:
        run_all_locus,
