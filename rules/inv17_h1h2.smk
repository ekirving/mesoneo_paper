#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand, temp, ancient

from scripts.utils import get_modern_samples, get_ancient_samples


rule bcftools_1000G:
    # 1000G phase 3 callset
    input:
        vcf=ancient(
            config["1000G"]["vcf_path"]
            + "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        ),
    output:
        vcf=temp("data/imputed_unfiltered/{dataset}-{population}-chr{chr}-1000G.bcf"),
        csi=temp("data/imputed_unfiltered/{dataset}-{population}-chr{chr}-1000G.bcf.csi"),
    params:
        samples=lambda wildcards: ",".join(get_modern_samples(config, wildcards)),
    threads: 4
    shell:
        "bcftools view"
        " --threads {threads}"
        " --samples {params.samples}"
        " --output-type b"
        " --output {output.vcf}"
        " {input.vcf} && "
        "bcftools index --threads {threads} {output.vcf}"


rule bcftools_unfiltered:
    # imputed best guess haplotypes (FORMAT/GT phased with | character) for all variants, no filtering
    input:
        vcf=config["imputed"]["vcf_path"] + "chr{chr}.haplotypes.unfiltered.N1664.D15062020.GLIMPSE.vcf.gz",
    output:
        vcf=temp("data/imputed_unfiltered/{dataset}-{population}-chr{chr}-unfiltered.bcf"),
        csi=temp("data/imputed_unfiltered/{dataset}-{population}-chr{chr}-unfiltered.bcf.csi"),
    params:
        samples=lambda wildcards: ",".join(get_ancient_samples(config, wildcards)),
    threads: 4
    shell:
        "bcftools view"
        " --threads {threads}"
        " --samples {params.samples}"
        " --output-type b"
        " --output {output.vcf}"
        " {input.vcf} && "
        "bcftools index --threads {threads} {output.vcf}"


rule bcftools_merge:
    input:
        vcf1="data/imputed_unfiltered/{dataset}-{population}-chr{chr}-1000G.bcf",
        csi1="data/imputed_unfiltered/{dataset}-{population}-chr{chr}-1000G.bcf.csi",
        vcf2="data/imputed_unfiltered/{dataset}-{population}-chr{chr}-unfiltered.bcf",
        csi2="data/imputed_unfiltered/{dataset}-{population}-chr{chr}-unfiltered.bcf.csi",
        vcf3=config["imputed"]["vcf_path"] + "chr{chr}.haplotypes.unfiltered.N1664.D15062020.GLIMPSE.vcf.gz",
    output:
        vcf=temp("data/imputed_unfiltered/{dataset}-{population}-chr{chr}-unfiltered-1000G.bcf"),
        csi=temp("data/imputed_unfiltered/{dataset}-{population}-chr{chr}-unfiltered-1000G.bcf.csi"),
    threads: 8
    shell:
        "bcftools merge"
        " --threads {threads}"
        " --regions-file {input.vcf3} "
        " --output-type b"
        " --output {output.vcf}"
        " {input.vcf1}"
        " {input.vcf2} && "
        "bcftools index --threads {threads} {output.vcf}"


rule bcftools_concat:
    input:
        expand(
            "data/imputed_unfiltered/{dataset}-{population}-chr{chr}-unfiltered-1000G.bcf", chr=config["chroms"], allow_missing=True
        ),
    output:
        vcf="data/imputed_unfiltered/{dataset}-{population}-chrALL-unfiltered-1000G.bcf",
        csi="data/imputed_unfiltered/{dataset}-{population}-chrALL-unfiltered-1000G.bcf.csi",
    threads: 88
    shell:
        "bcftools concat"
        " --threads {threads}"
        " --output-type b"
        " --output {output.vcf}"
        " {input} && "
        "bcftools index --threads {threads} {output.vcf}"


rule andres_report:
    input:
        data=config["andres"]["targets"],
        info="variants/{dataset}-{population}-info_info.tsv.gz",
    output:
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-andres_report.tsv",
    log:
        "clues/{dataset}-{population}-{mode}-{ancestry}-andres_report.tsv.log",
    shell:
        "python scripts/clues_report.py"
        " --data {input.data}"
        " --columns rsid"
        " --info {input.info}"
        " --dataset {wildcards.dataset}"
        " --population {wildcards.population}"
        " --mode {wildcards.mode}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.tsv} &> {log}"


rule inv17_h1h2_report:
    input:
        data=config["inv17_h1h2"]["targets"],
        info="variants/{dataset}-{population}-info_info.tsv.gz",
    output:
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-inv17_h1h2_report.tsv",
    log:
        "clues/{dataset}-{population}-{mode}-{ancestry}-inv17_h1h2_report.tsv.log",
    shell:
        "python scripts/clues_report.py"
        " --data {input.data}"
        " --columns SNP"
        " --info {input.info}"
        " --dataset {wildcards.dataset}"
        " --population {wildcards.population}"
        " --mode {wildcards.mode}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.tsv} &> {log}"
