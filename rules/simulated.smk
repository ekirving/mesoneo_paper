#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from scripts.utils import get_modern_samples, get_ancient_samples


rule simulated_ancient_mac_zero:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
    output:
        "variants/{dataset}-{population}-ancient_mac_zero.txt",
    params:
        samples=lambda wildcards: ",".join(get_ancient_samples(config, wildcards)),
    shell:
        r"bcftools view --samples {params.samples} --output-type u {input.vcf} | "
        r"bcftools view --max-ac 0 --output-type u | "
        r"bcftools query --format '%CHROM:%POS:%REF:%ALT\n' > {output}"


rule simulated_modern_freqs:
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
        zero="variants/{dataset}-{population}-ancient_mac_zero.txt",
    output:
        "variants/{dataset}-{population}-modern_freqs.tsv",
    params:
        samples=lambda wildcards: ",".join(get_modern_samples(config, wildcards)),
        columns=r"\t".join(["variant", "chr", "rsid", "ref", "alt", "AC", "AN"]),
    shell:
        r"printf '{params.columns}\n' > {output} && "
        r"bcftools view --samples {params.samples} --output-type u {input.vcf} | "
        r"bcftools +fill-tags --output-type u -- --tags AC,AN | "
        r"bcftools view --types snps --max-alleles 2 --min-ac 1 --output-type u | "
        r"bcftools annotate --set-id +'chr%CHROM:%POS' --output-type u | "
        r"bcftools query --format '%CHROM:%POS:%REF:%ALT\t%CHROM\t%ID\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\n' | "
        r"grep -vwF -f {input.zero} >> {output}"


checkpoint simulated_pair_snps:
    input:
        gwas="variants/ancestral_paths_v3-all-gwas.tsv",
        neut="variants/{dataset}-{population}-modern_freqs.tsv",
        anc="1000G/1000G_chrAll_ancestral.tsv.gz",
    output:
        protected("variants/{dataset}-{population}-pairs.tsv"),
    log:
        "variants/{dataset}-{population}-pairs.log",
    wildcard_constraints:
        dataset="|".join(["chr3_true_paths", "chr3_inferred_paths", "simulated_relate_painted"]),
    shell:
        "python scripts/neutral_pair_snps.py"
        " --gwas {input.gwas}"
        " --neutrals {input.neut}"
        " --ancestral {input.anc}"
        " --simulated "
        " --output {output} &> {log}"


ruleorder: simulated_pair_snps > neutral_pair_snps
