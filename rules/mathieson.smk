#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import expand

ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]


rule mathieson_manhattan_harvester:
    input:
        tsv=config["mathieson"]["targets"],
    output:
        tsv="mathieson/mathieson-harvester.tsv",
    log:
        "mathieson/mathieson-harvester.tsv.log",
    shell:
        "{config[harvester][path]}/harvester"
        " -chrcolumn 2"
        " -lcolumn 3"
        " -pcolumn 6"
        " -file {input.tsv}"
        " -out {output.tsv} &> {log}"


checkpoint mathieson_significant:
    input:
        tsv=config["mathieson"]["targets"],
    output:
        tsv="mathieson/mathieson-significant.tsv",
    shell:
        r"awk -F'\t' 'NR==1 || $6<5e-8 {{print $0}}' {input} > {output}"


rule mathieson_report_mode_ancestry:
    input:
        data="mathieson/mathieson-significant.tsv",
        info="variants/{dataset}-{population}-info_info.tsv.gz",
    output:
        tsv="mathieson/{dataset}-{population}-{mode}-{ancestry}-mathieson_report.tsv",
    log:
        "mathieson/{dataset}-{population}-{mode}-{ancestry}-mathieson_report.tsv.log",
    shell:
        "python scripts/clues_report.py"
        " --data {input.data}"
        " --columns ID"
        " --info {input.info}"
        " --dataset {wildcards.dataset}"
        " --population {wildcards.population}"
        " --mode {wildcards.mode}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.tsv} &> {log}"


rule mathieson_report:
    input:
        "mathieson/{dataset}-{population}-modern-ALL-mathieson_report.tsv",
        expand(
            "mathieson/{dataset}-{population}-ancient-{ancestry}-mathieson_report.tsv",
            ancestry=ANCESTRIES,
            allow_missing=True,
        ),
    output:
        "mathieson/{dataset}-{population}-mathieson_report.tsv",
    shell:
        "head -n1 {input[0]} > {output} && tail --quiet -n+2 {input} >> {output}"
