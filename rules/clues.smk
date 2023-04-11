#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

"""
CLUES analyses

Rules for inference of allele frequency trajectories and selection coefficients with CLUES.

See https://github.com/35ajstern/clues/
"""

import gzip
import json

from snakemake.io import expand, protected, temp, unpack, ancient

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
        coal=lambda wildcards: ancient(
            "relate/{panel}-{pops}-popsize.coal".format(
                panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
            )
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
        "python bin/clues/inference.py"
        " --lik"
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
        "python bin/clues/inference.py"
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
        "freq": ancient(
            expand("clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-ALL-any.freq", **params)
        ),
        "bins": expand("clues/{dataset}-{population}-time.bins", **params),
        "coal": expand("relate/{panel}-{pops}-popsize.coal", **params),
        "timeb": expand("relate/{panel}/{pops}/{rsid}/{panel}-{pops}-popsize-allsnps-{rsid}.timeb", **params),
        "meta": ancient(expand("variants/metadata/{reference}/{{rsid}}.json", **params)),
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
        "python bin/clues/inference.py"
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
    resources:
        matplotlib=1, # running too many `matplotlib` jobs simultaneously causes them to crash
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
        pairs="variants/{dataset}-{population}-pairs.tsv",
    output:
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-filtered-clues_report.tsv",
    wildcard_constraints:
        mode="ancient|modern"
    shell:
        "Rscript scripts/filter_refbias.R"
        " --data {input.tsv}"
        " --gwas {input.gwas}"
        " --neut {input.neut}"
        " --post {input.post}"
        " --pairs {input.pairs}"
        " --output {output.tsv}"


rule clues_simulation_report_mode_ancestry:
    input:
        pairs="variants/{dataset}-{population}-pairs.tsv",
    output:
        tsv="clues/{dataset}-{population}-simulated-{ancestry}-filtered-clues_report.tsv",
    log:
        "clues/{dataset}-{population}-simulated-{ancestry}-filtered-clues_report.tsv.log",
    wildcard_constraints:
        dataset="|".join(["chr3_true_paths", "chr3_inferred_paths", "simulated_relate_painted"]),
    shell:
        "python scripts/clues_report.py"
        " --data {input.pairs}"
        " --columns neutral"
        " --dataset {wildcards.dataset}"
        " --population {wildcards.population}"
        " --mode ancient"
        " --ancestry {wildcards.ancestry}"
        " --output {output.tsv} &> {log}"


rule clues_sweep_detection:
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
        tsv="clues/{dataset}-{population}-{mode}-{ancestry}-filtered-sweeps.tsv",
    shell:
        "Rscript scripts/clues_sweep_detection.R"
        " --data {input.tsv}"
        " --pairs {input.pairs}"
        " --unmapped {input.unmapped}"
        " --flipped {input.flipped}"
        " --output {output.tsv}"


rule clues_sweep_detection_merged:
    input:
        peaks_all="clues/{dataset}-{population}-{mode}-ALL-filtered-sweeps.tsv",
        peaks_whg="clues/{dataset}-{population}-{mode}-WHG-filtered-sweeps.tsv",
        peaks_ehg="clues/{dataset}-{population}-{mode}-EHG-filtered-sweeps.tsv",
        peaks_chg="clues/{dataset}-{population}-{mode}-CHG-filtered-sweeps.tsv",
        peaks_ana="clues/{dataset}-{population}-{mode}-ANA-filtered-sweeps.tsv",
    output:
        tsv="clues/{dataset}-{population}-{mode}-filtered-sweeps-merged.tsv",
    shell:
        "Rscript scripts/clues_sweep_detection_merged.R"
        " --peaks-all {input.peaks_all}"
        " --peaks-whg {input.peaks_whg}"
        " --peaks-ehg {input.peaks_ehg}"
        " --peaks-chg {input.peaks_chg}"
        " --peaks-ana {input.peaks_ana}"
        " --output {output.tsv}"


rule clues_report_filter_refbias:
    input:
        "clues/{dataset}-{population}-modern-ALL-filtered-clues_report.tsv",
        expand(
            "clues/{dataset}-{population}-ancient-{ancestry}-filtered-clues_report.tsv",
            ancestry=ANCESTRIES,
            allow_missing=True,
        ),
    output:
        "clues/{dataset}-{population}-filtered-clues_report.tsv",
    shell:
        "head -n1 {input[0]} > {output} && tail -q -n+2 {input} >> {output}"


rule clues_simulation_report:
    input:
        expand(
            "clues/{dataset}-{population}-simulated-{ancestry}-filtered-clues_report.tsv",
            ancestry=ANCESTRIES,
            allow_missing=True,
        ),
    output:
        "clues/{dataset}-{population}-simulated-clues_report.tsv",
    wildcard_constraints:
        dataset="|".join(["chr3_true_paths", "chr3_inferred_paths", "simulated_relate_painted"]),
    shell:
        "head -n1 {input[0]} > {output} && tail -q -n+2 {input} >> {output}"


rule clues_report_nearest_gene:
    input:
        tsv="clues/{dataset}-{population}-filtered-clues_report.tsv",
        bed=lambda wildcards: "bed/{reference}-gene_names.bed".format(
            reference=config["samples"][wildcards.dataset]["reference"]
        ),
    output:
        snps=temp("clues/{dataset}-{population}-filtered-clues_report.bed"),
        genes="clues/{dataset}-{population}-filtered-clues_report-genes.bed",
    shell:
        r"""awk -F'\t' 'NR>1 {{ print $3"\t"$4-1"\t"$5"\t"$1 }}' {input.tsv} | bedtools sort | uniq > {output.snps} && """
        r"""bedtools closest -a {output.snps} -b {input.bed} > {output.genes}"""


rule clues_plot_main_text_figure:
    input:
        data="clues/{dataset}-{population}-filtered-clues_report.tsv",
        genes="clues/{dataset}-{population}-filtered-clues_report-genes.bed",
        pairs="variants/{dataset}-{population}-pairs.tsv",
        unmapped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_unmapped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        flipped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_flipped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        peaks_all="clues/{dataset}-{population}-ancient-ALL-filtered-sweeps.tsv",
        peaks_whg="clues/{dataset}-{population}-ancient-WHG-filtered-sweeps.tsv",
        peaks_ehg="clues/{dataset}-{population}-ancient-EHG-filtered-sweeps.tsv",
        peaks_chg="clues/{dataset}-{population}-ancient-CHG-filtered-sweeps.tsv",
        peaks_ana="clues/{dataset}-{population}-ancient-ANA-filtered-sweeps.tsv",
        merged="clues/{dataset}-{population}-ancient-filtered-sweeps-merged.tsv",
        mathieson="mathieson/mathieson-sweeps.tsv",
    output:
        png="figs/{dataset}-{population}-filtered-main_figure.png",
    shell:
        "Rscript scripts/clues_plot_main_text_figure.R"
        " --data {input.data}"
        " --genes {input.genes}"
        " --pairs {input.pairs}"
        " --unmapped {input.unmapped}"
        " --flipped {input.flipped}"
        " --peaks-all {input.peaks_all}"
        " --peaks-whg {input.peaks_whg}"
        " --peaks-ehg {input.peaks_ehg}"
        " --peaks-chg {input.peaks_chg}"
        " --peaks-ana {input.peaks_ana}"
        " --merged {input.merged}"
        " --mathieson {input.mathieson}"
        " --output {output.png}"


rule clues_plot_extended_data_figures:
    input:
        data="clues/{dataset}-{population}-filtered-clues_report.tsv",
        genes="clues/{dataset}-{population}-filtered-clues_report-genes.bed",
        pairs="variants/{dataset}-{population}-pairs.tsv",
        unmapped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_unmapped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        flipped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_flipped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        peaks_all="clues/{dataset}-{population}-ancient-ALL-filtered-sweeps.tsv",
        peaks_whg="clues/{dataset}-{population}-ancient-WHG-filtered-sweeps.tsv",
        peaks_ehg="clues/{dataset}-{population}-ancient-EHG-filtered-sweeps.tsv",
        peaks_chg="clues/{dataset}-{population}-ancient-CHG-filtered-sweeps.tsv",
        peaks_ana="clues/{dataset}-{population}-ancient-ANA-filtered-sweeps.tsv",
        merged="clues/{dataset}-{population}-ancient-filtered-sweeps-merged.tsv",
        mathieson="mathieson/mathieson-sweeps.tsv",
    output:
        "figs/extended/{dataset}-{population}.done",
    params:
        prefix="figs/extended/{dataset}-{population}-locus",
    shell:
        "Rscript scripts/clues_plot_extended_data_figures.R"
        " --data {input.data}"
        " --genes {input.genes}"
        " --pairs {input.pairs}"
        " --unmapped {input.unmapped}"
        " --flipped {input.flipped}"
        " --peaks-all {input.peaks_all}"
        " --peaks-whg {input.peaks_whg}"
        " --peaks-ehg {input.peaks_ehg}"
        " --peaks-chg {input.peaks_chg}"
        " --peaks-ana {input.peaks_ana}"
        " --merged {input.merged}"
        " --mathieson {input.mathieson}"
        " --output {params.prefix} && "
        "touch {output}"


rule clues_plot_supplemental_figures:
    input:
        data="clues/{dataset}-{population}-filtered-clues_report.tsv",
        pairs="variants/{dataset}-{population}-pairs.tsv",
        unmapped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_unmapped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        flipped=lambda wildcards: "relate/{panel}-{pops}-popsize-allsnps_flipped.tsv.gz".format(
            panel="1000G_phase3", pops="_".join(get_modern_pops(config, wildcards))
        ),
        peaks_mod="clues/{dataset}-{population}-modern-ALL-filtered-sweeps.tsv",
        peaks_all="clues/{dataset}-{population}-ancient-ALL-filtered-sweeps.tsv",
        peaks_whg="clues/{dataset}-{population}-ancient-WHG-filtered-sweeps.tsv",
        peaks_ehg="clues/{dataset}-{population}-ancient-EHG-filtered-sweeps.tsv",
        peaks_chg="clues/{dataset}-{population}-ancient-CHG-filtered-sweeps.tsv",
        peaks_ana="clues/{dataset}-{population}-ancient-ANA-filtered-sweeps.tsv",
        mathieson="mathieson/mathieson-sweeps.tsv",
    output:
        "figs/supplement/{dataset}-{population}-modern.png",
        "figs/supplement/{dataset}-{population}-ancient.png",
        "figs/supplement/{dataset}-{population}-ancestral-gwas.png",
        "figs/supplement/{dataset}-{population}-ancestral-control.png",
    params:
        prefix="figs/supplement/{dataset}-{population}",
    shell:
        "Rscript scripts/clues_plot_supplement_figures.R"
        " --data {input.data}"
        " --pairs {input.pairs}"
        " --unmapped {input.unmapped}"
        " --flipped {input.flipped}"
        " --peaks-mod {input.peaks_mod}"
        " --peaks-all {input.peaks_all}"
        " --peaks-whg {input.peaks_whg}"
        " --peaks-ehg {input.peaks_ehg}"
        " --peaks-chg {input.peaks_chg}"
        " --peaks-ana {input.peaks_ana}"
        " --mathieson {input.mathieson}"
        " --output {params.prefix}"


rule clues_plot_simulated_figure:
    input:
        data="clues/{dataset}-{population}-simulated-clues_report.tsv",
        pairs="variants/{dataset}-{population}-pairs.tsv",
        peaks_all="clues/{dataset}-{population}-simulated-ALL-filtered-sweeps.tsv",
        peaks_whg="clues/{dataset}-{population}-simulated-WHG-filtered-sweeps.tsv",
        peaks_ehg="clues/{dataset}-{population}-simulated-EHG-filtered-sweeps.tsv",
        peaks_chg="clues/{dataset}-{population}-simulated-CHG-filtered-sweeps.tsv",
        peaks_ana="clues/{dataset}-{population}-simulated-ANA-filtered-sweeps.tsv",
        merged="clues/{dataset}-{population}-simulated-filtered-sweeps-merged.tsv",
    output:
        png="figs/supplement/{dataset}-{population}-simulated.png",
    wildcard_constraints:
        dataset="|".join(["chr3_true_paths", "chr3_inferred_paths", "simulated_relate_painted"]),
    shell:
        "Rscript scripts/clues_plot_simulated_figure.R "
        " --data {input.data}"
        " --pairs {input.pairs}"
        " --peaks-all {input.peaks_all}"
        " --peaks-whg {input.peaks_whg}"
        " --peaks-ehg {input.peaks_ehg}"
        " --peaks-chg {input.peaks_chg}"
        " --peaks-ana {input.peaks_ana}"
        " --merged {input.merged}"
        " --output {output.png}"
