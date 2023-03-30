#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"


rule plot_samples_map:
    input:
        tsv=lambda wildcards: config["samples"][wildcards.dataset]["metadata"],
    output:
        png="figs/{dataset}-samples_map.png",
    conda:
        "../envs/plot_samples_map.yaml"
    shell:
        "Rscript scripts/map_samples.R"
        " --samples {input.tsv}"
        " --output {output.png}"
