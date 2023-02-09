#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand, protected, temp, ancient


rule tgp_1000g_ancestral_chr:
    input:
        ancient(
            config["1000G"]["vcf_path"]
            + "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
        ),
    output:
        temp("1000G/1000G_chr{chr}_ancestral.tsv"),
    shell:
        r"printf 'variant\tancestral\n' > {output} && "
        r"bcftools view --types snps --max-alleles 2 --output-type u {input} | "
        r"bcftools query --format '%CHROM:%POS:%REF:%ALT\t%INFO/AA\n' | "
        r"sed 's/|||//' | tr [a-z] [A-Z] | grep -P '[ACGT]$' >> {output}"


rule tgp_1000g_ancestral:
    input:
        expand("1000G/1000G_chr{chr}_ancestral.tsv", chr=config["chroms"]),
    output:
        protected("1000G/1000G_chrAll_ancestral.tsv.gz"),
    shell:
        "cat {input} | gzip -c > {output}"
