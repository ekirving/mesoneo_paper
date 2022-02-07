#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

# pad all non-neutral loci by +/- 100 kb
NON_NEUTRAL_BUFFER = 50000


rule bedtools_genes:
    input:
        "ensembl/{reference}-annotation.gtf",
    output:
        "bed/{reference}-genes.bed",
    shell:
        'awk \'$3 == "gene" {{ print $1"\\t"$4-1"\\t"$5 }}\' {input} | '
        "bedtools sort | "
        "bedtools merge > {output}"


rule bedtools_variants:
    input:
        "variants/{reference}_report.tsv",
    output:
        "bed/{reference}-variants.bed",
    shell:
        'awk \'$2 == "SNP" {{ print $3"\\t"$4-1"\\t"$5}}\' {input} | '
        "bedtools sort | "
        "bedtools merge > {output}"


rule bedtools_neutral_loci:
    input:
        genes="bed/{reference}-genes.bed",
        snps="bed/{reference}-variants.bed",
        genome="data/{reference}.size",
    output:
        "bed/{reference}-neutral.bed",
    params:
        buffer=NON_NEUTRAL_BUFFER,
    shell:
        "cat {input.genes} {input.snps} | "
        "grep -P '^(\d+|X|Y)' | "
        "bedtools slop -b {params.buffer} -g {input.genome} | "
        "bedtools sort -g {input.genome} | "
        "bedtools merge | "
        "bedtools complement -i stdin -g {input.genome} > {output}"


rule bedtools_neutral_strict:
    input:
        bed="bed/{reference}-neutral.bed",
        mask=config["1000G"]["mask_path"] + "20141020.strict_mask.whole_genome.bed",
    output:
        "bed/{reference}-neutral_strict.bed",
    shell:
        "bedtools intersect -a {input.bed} -b <(sed 's/chr//' {input.mask}) > {output}"
