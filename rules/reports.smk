#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import ancient


rule bcftools_snp_freqs:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
        bed=lambda wildcards: "bed/{reference}-variants.bed".format(
            reference=config["samples"][wildcards.dataset]["reference"]
        ),
    output:
        "variants/{dataset}-{population}-freqs.tsv",
    params:
        columns=r"\t".join(["chr", "pos", "rsid", "ref", "alt", "AC", "AN"]),
    shell:
        r"printf '{params.columns}\n' > {output} && "
        r"bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n'"
        r" --regions-file {input.bed} {input.vcf} >> {output}"


checkpoint bcftools_info_info:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
    output:
        "variants/{dataset}-{population}-info_info.tsv.gz",
    shell:
        r"printf 'rsid\tchrom\tstart\tinfo\n'| gzip -c > {output} && "
        r"bcftools query --format '%ID\t%CHROM\t%POS\t%INFO/INFO\n' {input} | gzip -c >> {output}"


rule variant_report_population:
    input:
        report=lambda wildcards: "variants/{reference}_report.tsv".format(
            reference=config["samples"][wildcards.dataset]["reference"]
        ),
        freqs="variants/{dataset}-{population}-freqs.tsv",
    output:
        bed="variants/{dataset}-{population}_report.bed",
    shell:
        "python scripts/variant_report_population.py"
        " --report {input.report}"
        " --freqs {input.freqs}"
        " --output {output}"


checkpoint variant_strict_report:
    input:
        bed="variants/{dataset}-{population}_report.bed",
        mask=ancient(config["1000G"]["mask_path"] + "20141020.strict_mask.whole_genome.bed"),
    output:
        bed="variants/{dataset}-{population}_report_strict.bed",
    shell:
        "bedtools intersect"
        " -a {input.bed}"
        " -b <(sed 's/chr//' {input.mask}) > {output.bed}"
