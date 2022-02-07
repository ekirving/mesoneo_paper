#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


rule ensembl_variation:
    # we use a customised version of the ensembl-variation wrapper to fetch the GRCh37 build of Ensembl
    # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/reference/ensembl-variation.html
    output:
        vcf="ensembl/{reference}-variation.vcf.gz",
    params:
        species="Homo_sapiens",
        release="99",
        type="all",
    log:
        "ensembl/{reference}-variation.log",
    script:
        "../scripts/ensembl_variation.py"


rule ensembl_variation_index:
    input:
        vcf="ensembl/{reference}-variation.vcf.gz",
    output:
        csi="ensembl/{reference}-variation.vcf.gz.csi",
    shell:
        "bcftools index {input}"


rule ensembl_annotation:
    # we use a customised version of the ensembl-annotation wrapper to fetch the GRCh37 build of Ensembl
    # https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/reference/ensembl-annotation.html
    output:
        "ensembl/{reference}-annotation.gtf",
    params:
        species="Homo_sapiens", # release 87 is the last version of the GRCh37 annotation
        release=lambda wildcards: "87" if wildcards.reference == "GRCh37" else "99",
        build=lambda wildcards: wildcards.reference,
        fmt="gtf",
    log:
        "ensembl/{reference}-annotation.log",
    script:
        "../scripts/ensembl_annotation.py"
