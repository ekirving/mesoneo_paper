#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import ancient

# maximum concurrent requests
ENSEMBL_MAX_CONCURRENT = 15


rule variant_ensembl:
    output:
        var="variants/ensembl/{reference}/var/{rsid}.json",
        vep="variants/ensembl/{reference}/vep/{rsid}.json",
    resources:
        ensembl_api=1,
    shell:
        "python scripts/variant_ensembl.py"
        " --ref {wildcards.reference}"
        " --rsid {wildcards.rsid}"
        " --var {output.var}"
        " --vep {output.vep}"


checkpoint variant_metadata:
    input:
        var="variants/ensembl/{reference}/var/{rsid}.json",
        vep="variants/ensembl/{reference}/vep/{rsid}.json",
        gwas=ancient("gwascat/gwas_catalog_targets.tsv"),
    output:
        "variants/metadata/{reference}/{rsid}.json",
    resources:
        mem_mb=8000,
    wildcard_constraints:
        rsid="rs\d+",
    shell:
        "python scripts/variant_metadata.py"
        " --rsid {wildcards.rsid}"
        " --var {input.var}"
        " --vep {input.vep}"
        " --gwas {input.gwas}"
        " --out {output}"


checkpoint variant_simulated:
    output:
        "variants/metadata/{reference}/chr{chr}:{pos}.json",
    params:
        data=lambda wildcards: {
            "rsid": "chr" + wildcards.chr + ":" + wildcards.pos,
            "type": "SNP",
            "chrom": wildcards.chr,
            "start": wildcards.pos,
            "end": wildcards.pos,
            "allele": "1",
            "alleles": "0/1",
            "ancestral": "0",
            "derived": "1",
        },
    run:
        import json

        with open(output[0], "w") as fout:
            json.dump(params.data, fout, indent=2)


rule variant_label:
    input:
        meta=lambda wildcards: "variants/metadata/{reference}/{{rsid}}.json".format(
            reference=config["samples"][wildcards.dataset]["reference"]
        ),
    output:
        "variants/labels/{dataset}/{population}/{dataset}-{population}-{rsid}-label.json",
    params:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
    shell:
        "python scripts/variant_label.py"
        " --vcf {params.vcf}"
        " --meta {input.meta}"
        " --out {output}"


rule variant_batch_all:
    input:
        gwas="gwascat/gwas_catalog_targets.tsv",
    output:
        log="variants/metadata/{reference}/batch.log",
    params:
        ref="{reference}",
        var="variants/ensembl/{reference}/var/[].json",
        vep="variants/ensembl/{reference}/vep/[].json",
        met="variants/metadata/{reference}/[].json",
    threads: ENSEMBL_MAX_CONCURRENT
    shell:
        "mkdir -p $(dirname {params.var}) && "
        "mkdir -p $(dirname {params.vep}) && "
        "(awk 'NR>1 {{print $1}}' {input.gwas} | sort -u > {output.log}) && "
        "(parallel -P {threads} -a {output.log} -I [] '"
        " if [[ ! -f {params.var} ]] || [[ ! -f {params.vep} ]]; then"
        "  python scripts/variant_ensembl.py --ref {params.ref} --rsid [] --var {params.var} --vep {params.vep}; fi && "
        " if [[ ! -f {params.met} ]]; then"
        "  python scripts/variant_metadata.py --rsid [] --var {params.var} --vep {params.vep} --gwas {input.gwas} --out {params.met}; fi;"
        "')"


rule variant_report_metadata:
    input:
        "variants/metadata/{reference}/batch.log",
    output:
        "variants/{reference}_report.tsv",
    shell:
        "python scripts/variant_report_metadata.py"
        " --ref {wildcards.reference}"
        " --rsids {input}"
        " --out {output}"
