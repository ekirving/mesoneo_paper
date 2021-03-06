#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import gzip
import json

from snakemake.io import expand, protected, temp, unpack

from scripts.utils import trim_ext, get_modern_pops

ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]


rule clues_ancient_samples:
    input:
        meta=lambda wildcards: "variants/metadata/{reference}/{{rsid}}.json".format(
            reference=config["samples"][wildcards.dataset]["reference"]
        ),
    output:
        anc="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{ancestry}-{sex}.ancient",
        freq="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{ancestry}-{sex}.freq",
    log:
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{ancestry}-{sex}.ancient.log",
    params:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
        meta=lambda wildcards, input: json.load(open(str(input.meta))),
        ancestral=lambda wildcards, input: json.load(open(str(input.meta))).get("ancestral") or "N",
    shell:
        "python scripts/clues_ancient_samples.py"
        " --vcf {params.vcf}"
        " --chr {params.meta[chrom]}"
        " --pos {params.meta[start]}"
        " --ancestral {params.ancestral}"
        " --dataset {wildcards.dataset}"
        " --population {wildcards.population}"
        " --ancestry {wildcards.ancestry}"
        " --sex {wildcards.sex}"
        " --gen-time {config[gen_time]}"
        " --mod-freq {output.freq}"
        " --output {output.anc} 2> {log}"


rule clues_time_bins:
    output:
        bins="clues/{dataset}-{population}-time.bins",
    shell:
        "python scripts/clues_time_bins.py"
        " --dataset {wildcards.dataset}"
        " --population {wildcards.population}"
        " --gen-time {config[gen_time]}"
        " --output {output}"


rule clues_inference_ancient:
    input:
        anc="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{ancestry}-{sex}.ancient",
        freq="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{ancestry}-{sex}.freq",
        bins="clues/{dataset}-{population}-time.bins",
        coal=lambda wildcards: "relate/{panel}-{pops}-popsize.coal".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
    output:
        epochs="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-{ancestry}-{sex}.epochs.npy",
        freqs="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-{ancestry}-{sex}.freqs.npy",
        post="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-{ancestry}-{sex}.post.npy",
    log:
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ancient-{ancestry}-{sex}.epochs.log",
    params:
        output=lambda wildcards, input, output: trim_ext(output.epochs, 2),
        daf=lambda wildcards, input: open(str(input.freq)).read().strip(),
        flag=lambda wildcards: "ancientSamps" if wildcards.ancestry == "ALL" else "ancientHaps",
    shell:
        "python {config[clues][path]}/inference.py"
        " --popFreq {params.daf}"
        " --coal {input.coal}"
        " --{params.flag} {input.anc}"
        " --timeBins {input.bins}"
        " --betaParam 0.5"
        " --out {params.output} &> {log}"


def clues_inference_modern_inputs(wildcards):
    """Get the Relate trees for this dataset"""
    inputs = clues_inference_both_inputs(wildcards)
    inputs.pop("anc")

    return inputs


def get_relate_allele(wildcards, input):
    """
    We need to flip the derived allele frequency if Relate has flipped the SNP during inference.
    """
    # get the list of SNPs that Relate flipped
    # noinspection PyUnresolvedReferences
    flipped_file = checkpoints.relate_is_flipped.get(panel="1000G_phase3", pops="FIN_GBR_TSI")

    with gzip.open(flipped_file.output[0], "r") as fin:
        flipped_snps = [line.decode("utf-8") for line in fin.read().split()]

    if wildcards.rsid in flipped_snps:
        # return the ancestral, so CLUES flips the trajectory
        return json.load(open(str(input.meta))).get("ancestral")
    else:
        return json.load(open(str(input.meta))).get("derived")

    return daf


rule clues_inference_modern:
    input:
        unpack(clues_inference_modern_inputs),
    output:
        epochs=protected(
            "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-modern-ALL-any.epochs.npy"
        ),
        freqs=protected("clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-modern-ALL-any.freqs.npy"),
        post=protected("clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-modern-ALL-any.post.npy"),
    log:
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-modern-ALL-any.epochs.log",
    params:
        input=lambda wildcards, input: trim_ext(input.timeb),
        output=lambda wildcards, input, output: trim_ext(output.epochs, 2),
        daf=lambda wildcards, input: open(str(input.freq)).read().strip(),
        allele=get_relate_allele,
    shell:
        "python {config[clues][path]}/inference.py"
        " --A1 {params.allele}"
        " --times {params.input}"
        " --coal {input.coal}"
        " --popFreq {params.daf}"
        " --timeBins {input.bins}"
        " --betaParam 0.5"
        " --out {params.output} &> {log}"


def clues_inference_both_inputs(wildcards):
    """Get the Relate trees for this dataset"""
    params = dict(**wildcards)
    params["reference"] = config["samples"][wildcards.dataset]["reference"]
    params["panel"] = "1000G_phase3"
    params["pops"] = "_".join(get_modern_pops(config, wildcards))

    return {
        "anc": expand("clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ALL-any.ancient", **params),
        "freq": expand("clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ALL-any.freq", **params),
        "bins": expand("clues/{dataset}-{population}-time.bins", **params),
        "coal": expand("relate/{panel}-{pops}-popsize.coal", **params),
        "timeb": expand("relate/{panel}/{pops}/{rsid}/{panel}-{pops}-popsize-allsnps-{rsid}.timeb", **params),
        "meta": expand("variants/metadata/{reference}/{{rsid}}.json", **params),
    }


rule clues_inference_both:
    input:
        unpack(clues_inference_both_inputs),
    output:
        epochs="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-both-ALL-any.epochs.npy",
        freqs="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-both-ALL-any.freqs.npy",
        post="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-both-ALL-any.post.npy",
    log:
        "clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-both-ALL-any.epochs.log",
    params:
        input=lambda wildcards, input: trim_ext(input.timeb),
        output=lambda wildcards, input, output: trim_ext(output.epochs, 2),
        daf=lambda wildcards, input: open(str(input.freq)).read().strip(),
        allele=get_relate_allele,
    shell:
        "python {config[clues][path]}/inference.py"
        " --A1 {params.allele}"
        " --times {params.input}"
        " --coal {input.coal}"
        " --popFreq {params.daf}"
        " --ancientSamps {input.anc}"
        " --timeBins {input.bins}"
        " --betaParam 0.5"
        " --out {params.output} &> {log}"


rule clues_parse_log:
    input:
        log="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{mode}-{ancestry}-{sex}.epochs.log",
    output:
        params="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{mode}-{ancestry}-{sex}.json",
    shell:
        "python scripts/clues_parse_log.py"
        " --rsid {wildcards.rsid}"
        " --mode {wildcards.mode}"
        " --ancestry {wildcards.ancestry}"
        " --sex {wildcards.sex}"
        " --log {input.log}"
        " --out {output.params}"


rule clues_plot_trajectory:
    input:
        epochs="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{mode}-{ancestry}-{sex}.epochs.npy",
        freqs="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{mode}-{ancestry}-{sex}.freqs.npy",
        post="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{mode}-{ancestry}-{sex}.post.npy",
        params="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{mode}-{ancestry}-{sex}.json",
        label="variants/labels/{dataset}/{population}/{dataset}-{population}-{rsid}-label.json",
    output:
        png="clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{mode}-{ancestry}-{sex}.png",
    params:
        input=lambda wildcards, input: trim_ext(input[0], 2),
        output=lambda wildcards, input, output: trim_ext(output[0]),
    shell:
        "python scripts/clues_plot_trajectory.py"
        " --gen-time {config[gen_time]}"
        " --params {input.params}"
        " --label {input.label}"
        " --ancestry {wildcards.ancestry}"
        " --sex {wildcards.sex}"
        " --ext png"
        " {params.input}"
        " {params.output}"


rule clues_report_mode_ancestry:
    input:
        pairs="variants/{dataset}-{population}-pairs.tsv",
        info="variants/{dataset}-{population}-info_info.tsv.gz",
    output:
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-clues_report.tsv",
    log:
        "clues/{dataset}-{population}-{mode}-{ancestry}-clues_report.tsv.log",
    shell:
        "python scripts/clues_report.py"
        " --data {input.pairs}"
        " --columns gwas,neutral"
        " --info {input.info}"
        " --dataset {wildcards.dataset}"
        " --population {wildcards.population}"
        " --mode {wildcards.mode}"
        " --ancestry {wildcards.ancestry}"
        " --output {output.tsv} &> {log}"


rule clues_report_mode_ancestry_filter_refbias:
    input:
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-clues_report.tsv",
        gwas="refbias/Evan_GWAS_Anc_1000g.txt.gz",
        neut="refbias/Evan_NEUTRAL_Anc_1000g.txt.gz",
        post="refbias/posterior-diff.tsv.gz",
    output:
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-filtered-clues_report.tsv",
    shell:
        "Rscript scripts/filter_refbias.R"
        " --data {input.tsv}"
        " --gwas {input.gwas}"
        " --neut {input.neut}"
        " --post {input.post}"
        " --output {output.tsv}"


rule clues_harvester_input:
    input:
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-filtered-clues_report.tsv",
        pairs="variants/{dataset}-{population}-pairs.tsv",
        unmapped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_unmapped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        flipped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_flipped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
    output:
        tsv=temp("clues/{dataset}-{population}-{mode}-{ancestry}-filtered-harvester-input.tsv"),
    shell:
        "Rscript scripts/harvester_input.R"
        " --data {input.tsv}"
        " --pairs {input.pairs}"
        " --unmapped {input.unmapped}"
        " --flipped {input.flipped}"
        " --output {output.tsv}"


rule clues_manhattan_harvester:
    input:
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-filtered-harvester-input.tsv",
    output:
        tmp=temp("clues/{dataset}-{population}-{mode}-{ancestry}-filtered-harvester-output.tsv"),
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-filtered-harvester.tsv",
    log:
        "clues/{dataset}-{population}-{mode}-{ancestry}-filtered-harvester.tsv.log",
    shell:
        "{config[harvester][path]}/harvester"
        " -chrcolumn 1"
        " -lcolumn 2"
        " -pcolumn 3"
        " -inlimit 0.01"
        " -peak-limit 5.87"
        " -file {input}"
        " -out {output.tmp} &> {log} && "
        "awk -F'\\t' 'BEGIN {{OFS=FS}} NR>1 {{split($3, range, \"-\"); $3=range[1]*10\"-\"range[2]*10}} {{print $0}}'"
        " {output.tmp} > {output.tsv}"


rule clues_report:
    input:
        "clues/{dataset}-{population}-modern-ALL-clues_report.tsv",
        expand(
            "clues/{dataset}-{population}-ancient-{ancestry}-clues_report.tsv", ancestry=ANCESTRIES, allow_missing=True
        ),
    output:
        "clues/{dataset}-{population}-clues_report.tsv",
    shell:
        "head -n1 {input[0]} > {output} && tail --quiet -n+2 {input} >> {output}"


rule clues_report_filter_refbias:
    input:
        expand(
            "clues/{dataset}-{population}-ancient-{ancestry}-filtered-clues_report.tsv",
            ancestry=ANCESTRIES,
            allow_missing=True,
        ),
    output:
        "clues/{dataset}-{population}-filtered-clues_report.tsv",
    shell:
        "head -n1 {input[0]} > {output} && tail --quiet -n+2 {input} >> {output}"


rule clues_pvalue:
    input:
        "clues/{dataset}-{population}-clues_report.tsv",
    output:
        "clues/{dataset}-{population}-clues_report_pvalue.tsv",
    shell:
        "Rscript scripts/clues_pvalue.R"
        " --data {input}"
        " --output {output}"


rule clues_plot_manhattan:
    input:
        data="clues/{dataset}-{population}-clues_report_pvalue.tsv",
        pairs="variants/{dataset}-{population}-pairs.tsv",
        known="data/mathieson/41586_2015_BFnature16152_MOESM270_ESM.txt",
        unmapped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_unmapped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        flipped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_flipped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
    output:
        pdf="clues/{dataset}-{population}-manhattan-{x}_vs_{y}.pdf",
    shell:
        "Rscript scripts/clues_plot_manhattan.R"
        " --data {input.data}"
        " --pairs {input.pairs}"
        " --known {input.known}"
        " --unmapped {input.unmapped}"
        " --flipped {input.flipped}"
        " --output {output.pdf}"
        " --facet-x {wildcards.x}"
        " --facet-y {wildcards.y}"


rule clues_plot_manhattan_chr3:
    input:
        data="clues/chr3_{paths}_paths-{population}-clues_report_pvalue.tsv",
        pairs="variants/chr3_simulated-pairs.tsv",
        known="data/mathieson/41586_2015_BFnature16152_MOESM270_ESM.txt",
    output:
        pdf="clues/chr3_{paths}_paths-{population}-manhattan-{x}_vs_{y}.pdf",
    shell:
        "Rscript scripts/clues_plot_manhattan.R"
        " --data {input.data}"
        " --pairs {input.pairs}"
        " --known {input.known}"
        " --output {output.pdf}"
        " --facet-x {wildcards.x}"
        " --facet-y {wildcards.y}"


ruleorder: clues_plot_manhattan_chr3 > clues_plot_manhattan


rule clues_plot_main_text_figure:
    input:
        data="clues/{dataset}-{population}-filtered-clues_report.tsv",
        pairs="variants/{dataset}-{population}-pairs.tsv",
        unmapped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_unmapped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        flipped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_flipped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        peaks_all="clues/{dataset}-{population}-ancient-ALL-filtered-harvester.tsv",
        peaks_whg="clues/{dataset}-{population}-ancient-WHG-filtered-harvester.tsv",
        peaks_ehg="clues/{dataset}-{population}-ancient-EHG-filtered-harvester.tsv",
        peaks_chg="clues/{dataset}-{population}-ancient-CHG-filtered-harvester.tsv",
        peaks_ana="clues/{dataset}-{population}-ancient-ANA-filtered-harvester.tsv",
        mathieson="mathieson/mathieson-harvester.tsv",
    output:
        pdf="figs/{dataset}-{population}-filtered-main_figure.png",
    shell:
        "Rscript scripts/clues_plot_main_text_figure.R"
        " --data {input.data}"
        " --pairs {input.pairs}"
        " --unmapped {input.unmapped}"
        " --flipped {input.flipped}"
        " --peaks-all {input.peaks_all}"
        " --peaks-whg {input.peaks_whg}"
        " --peaks-ehg {input.peaks_ehg}"
        " --peaks-chg {input.peaks_chg}"
        " --peaks-ana {input.peaks_ana}"
        " --mathieson {input.mathieson}"
        " --output {output.pdf}"
