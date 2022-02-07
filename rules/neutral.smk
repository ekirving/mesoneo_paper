#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import protected

from scripts.utils import get_modern_samples, get_modern_pops

# pad all non-neutral loci by +/- 50 kb
NON_NEUTRAL_BUFFER = 50000


rule bcftools_gwas_snps:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
        bed=lambda wildcards: "bed/{reference}-variants.bed".format(
            reference=config["samples"][wildcards.dataset]["reference"]
        ),
    output:
        "variants/{dataset}-{population}-gwas.tsv",
    params:
        samples=lambda wildcards: ",".join(get_modern_samples(config, wildcards)),
        columns=r"\t".join(["variant", "chr", "rsid", "ref", "alt", "AC", "AN"]),
    shell:
        r"printf '{params.columns}\n' > {output} && "
        r"bcftools view --types snps --max-alleles 2 --regions-file {input.bed} --samples {params.samples}"
        r" --output-type u {input.vcf} | "
        r"bcftools query --format '%CHROM:%POS:%REF:%ALT\t%CHROM\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' >> {output}"


rule bcftools_neutral_snps:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
        bed=lambda wildcards: "bed/{reference}-neutral_strict.bed".format(
            reference=config["samples"][wildcards.dataset]["reference"]
        ),
    output:
        "variants/{dataset}-{population}-neutrals.tsv",
    params:
        samples=lambda wildcards: ",".join(get_modern_samples(config, wildcards)),
        columns=r"\t".join(["variant", "chr", "rsid", "ref", "alt", "AC", "AN"]),
    shell:
        r"printf '{params.columns}\n' > {output} && "
        r"bcftools view --types snps --max-alleles 2 --regions-file {input.bed} --samples {params.samples}"
        r" --output-type u {input.vcf} | "
        r"bcftools query --format '%CHROM:%POS:%REF:%ALT\t%CHROM\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' >> {output}"


checkpoint neutral_pair_snps:
    input:
        gwas="variants/{dataset}-{population}-gwas.tsv",
        neut="variants/{dataset}-{population}-neutrals.tsv",
        anc="1000G/1000G_chrAll_ancestral.tsv.gz",
        bad=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_unmapped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
    output:
        protected("variants/{dataset}-{population}-pairs.tsv"),
    log:
        "variants/{dataset}-{population}-pairs.log",
    shell:
        "python scripts/neutral_pair_snps.py"
        " --gwas {input.gwas}"
        " --neutrals {input.neut}"
        " --ancestral {input.anc}"
        " --unmapped {input.bad}"
        " --output {output} &> {log}"
