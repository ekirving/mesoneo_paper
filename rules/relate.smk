#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json
from multiprocessing import cpu_count

from psutil import virtual_memory
from snakemake.io import directory, expand, protected, temp, unpack, ancient

from scripts.utils import trim_ext

# number of CPU cores
MAX_THREADS = cpu_count()

# find the number of threads nearest to 20 that divide evenly into the total cores
RELATE_NUM_THREADS = int(MAX_THREADS / max(1, round(MAX_THREADS / 20)))

# maximum of 90% of the RAM on current machine
MAX_MEM_MB = (virtual_memory().total / (1024 ** 2)) * 0.9

# memory usage is per thread, so give each thread the maximum amount of RAM
RELATE_MEM_MB = int((MAX_MEM_MB / MAX_THREADS) * RELATE_NUM_THREADS)

RELATE_FORMAT_ANCMUT = "a"
RELATE_FORMAT_NEWICK = "n"
RELATE_FORMAT_BINARY = "b"

RELATE_NO_FIXED_SNPS = 0
RELATE_ALL_SNPS = 1


rule relate_convert_vcf:
    input:
        vcf=(
            config["1000G"]["vcf_path"] + "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz"
        ),
    output:
        haps=temp("haps/{panel}_chr{chr}_convert.haps"),
        smpl="haps/{panel}_chr{chr}.sample",
    log:
        "haps/{panel}_chr{chr}_convert.log",
    params:
        stub=lambda wildcards, input: input.vcf.replace(".vcf.gz", ""),
    shell:
        "{config[relate][path]}/bin/RelateFileFormats"
        " --mode ConvertFromVcf"
        " --haps {output.haps}"
        " --sample {output.smpl}"
        " --input {params.stub} &> {log}"


rule relate_remove_non_biallelic_snps:
    input:
        haps="haps/{panel}_chr{chr}_convert.haps",
        smpl="haps/{panel}_chr{chr}.sample",
    output:
        haps=temp("haps/{panel}_chr{chr}_biallelic.haps"),
    log:
        "haps/{panel}_chr{chr}_biallelic.log",
    params:
        stub=lambda wildcards, input, output: trim_ext(output[0]),
    shell:
        "{config[relate][path]}/bin/RelateFileFormats"
        " --mode RemoveNonBiallelicSNPs"
        " --haps {input.haps}"
        " --output {params.stub} &> {log}"


rule relate_flip_haps_using_ancestor:
    input:
        haps="haps/{panel}_chr{chr}_biallelic.haps",
        smpl="haps/{panel}_chr{chr}.sample",
        anc=ancient(config["1000G"]["anc_path"] + "{chr}.humanc_e71.fa"),
    output:
        haps=temp("haps/{panel}_chr{chr}_ancestral.haps"),
    log:
        "haps/{panel}_chr{chr}_ancestral.log",
    params:
        stub=lambda wildcards, input, output: trim_ext(output[0]),
    shell:
        "{config[relate][path]}/bin/RelateFileFormats"
        " --mode FlipHapsUsingAncestor"
        " --haps {input.haps}"
        " --sample {input.smpl}"
        " --ancestor {input.anc}"
        " --output {params.stub} &> {log}"


rule relate_filter_haps_using_mask:
    input:
        haps="haps/{panel}_chr{chr}_ancestral.haps",
        smpl="haps/{panel}_chr{chr}.sample",
        mask=config["1000G"]["mask_path"] + "StrictMask/20140520.chr{chr}.strict_mask.fasta",
    output:
        haps="haps/{panel}_chr{chr}.haps",
        dist="haps/{panel}_chr{chr}.dist",
    log:
        "haps/{panel}_chr{chr}_mask.log",
    params:
        stub=lambda wildcards, input, output: trim_ext(output[0]),
    shell:
        "{config[relate][path]}/bin/RelateFileFormats"
        " --mode FilterHapsUsingMask"
        " --haps {input.haps}"
        " --sample {input.smpl}"
        " --mask {input.mask}"
        " --output {params.stub} &> {log}"


rule relate_generate_snp_annotations:
    input:
        haps="haps/{panel}_chr{chr}.haps",
        smpl="haps/{panel}_chr{chr}.sample",
        anc=ancient(config["1000G"]["anc_path"] + "{chr}.humanc_e71.fa"),
    output:
        annot="haps/{panel}_chr{chr}.annot",
    log:
        "haps/{panel}_chr{chr}_annot.log",
    params:
        stub=lambda wildcards, input, output: trim_ext(output[0]),
    shell:
        "{config[relate][path]}/bin/RelateFileFormats"
        " --mode GenerateSNPAnnotations"
        " --haps {input.haps}"
        " --sample {input.smpl}"
        " --poplabels {config[1000G][pop_labels]}"
        " --ancestor {input.anc}"
        " --output {params.stub} &> {log}"


rule relate_remove_samples:
    input:
        haps="haps/{panel}_chr{chr}.haps",
        smpl="haps/{panel}_chr{chr}.sample",
    output:
        haps="haps/{panel}-{pops}_chr{chr}.haps",
        smpl="haps/{panel}-{pops}_chr{chr}.sample",
        remv="haps/{panel}-{pops}_chr{chr}.remove",
        pops="haps/{panel}-{pops}_chr{chr}.poplabels",
    log:
        "haps/{panel}-{pops}_chr{chr}.log",
    params:
        stub=lambda wildcards, input, output: trim_ext(output[0]),
        pops=lambda wildcards: wildcards.pops.replace("_", "|"),
    shell:
        "grep -vwP '{params.pops}' {config[1000G][pop_labels]} | awk 'NR>1 {{print $1}}' > {output.remv} && "
        "{config[relate][path]}/bin/RelateFileFormats"
        " --mode RemoveSamples"
        " --haps {input.haps}"
        " --sample {input.smpl}"
        " --poplabels {config[1000G][pop_labels]}"
        " --flag {RELATE_ALL_SNPS}"
        " --input {output.remv}"
        " --output {params.stub} &> {log}"


rule relate_parallel:
    input:
        haps="haps/{panel}_chr{chr}.haps",
        smpl="haps/{panel}_chr{chr}.sample",
        map=config["1000G"]["map_path"] + "genetic_map_chr{chr}_combined_b37.txt",
        dist="haps/{panel}_chr{chr}.dist",
        annot="haps/{panel}_chr{chr}.annot",
    output:
        anc=protected("relate/{panel}_chr{chr}.anc.gz"),
        mut=protected("relate/{panel}_chr{chr}.mut.gz"),
        tmp=temp(directory(".tmp_{panel}_chr{chr}")),
    log:
        "relate/{panel}_chr{chr}.log",
    resources:
        mem_mb=RELATE_MEM_MB,
    threads: RELATE_NUM_THREADS
    params:
        mem_gb=lambda wildcards, input, output, resources, threads: int(resources.mem_mb / 1024 / threads), # mem_gb is per thread
        stub="{panel}_chr{chr}",
    shell:
        "mkdir -p {output.tmp} && cd {output.tmp} && "
        "{config[relate][path]}/scripts/RelateParallel/RelateParallel.sh"
        " --haps ../{input.haps}"
        " --sample ../{input.smpl}"
        " --map {input.map}"
        " --mu {config[relate][mu]}"
        " --Ne {config[relate][Ne]}"
        " --dist ../{input.dist}"
        " --annot ../{input.annot}"
        " --memory {params.mem_gb}"
        " --seed {wildcards.chr}"
        " --threads {threads}"
        " --output {params.stub} &> ../{log} && "
        "gzip {params.stub}.anc && "
        "gzip {params.stub}.mut && "
        "mv {params.stub}.* ../relate/"


rule relate_subtrees_for_subpopulation:
    input:
        anc="relate/{panel}_chr{chr}.anc.gz",
        mut="relate/{panel}_chr{chr}.mut.gz",
    output:
        anc="relate/{panel}-{pops}_chr{chr}.anc.gz",
        mut="relate/{panel}-{pops}_chr{chr}.mut.gz",
        lbl="relate/{panel}-{pops}_chr{chr}.poplabels",
    log:
        "relate/{panel}-{pops}_chr{chr}.log",
    params:
        pops=lambda wildcards: wildcards.pops.replace("_", ","),
        stub=lambda wildcards, input, output: trim_ext(output[0], 2),
    shell:
        "{config[relate][path]}/bin/RelateExtract"
        " --mode SubTreesForSubpopulation"
        " --anc {input.anc}"
        " --mut {input.mut}"
        " --poplabels {config[1000G][pop_labels]}"
        " --pop_of_interest {params.pops}"
        " --output {params.stub} &> {log} && "
        "gzip {params.stub}.anc && "
        "gzip {params.stub}.mut"


rule relate_popfile_for_metapopulation:
    input:
        "relate/{panel}-{pops}_chr1.poplabels",
    output:
        "relate/{panel}-{pops}.poplabels",
    shell:
        "awk 'NR>1 {{$2=$3}} $0' {input} > {output}"


rule relate_estimate_population_size:
    input:
        expand(
            ["relate/{panel}-{pops}_chr{chr}.anc.gz", "relate/{panel}-{pops}_chr{chr}.mut.gz"],
            chr=config["chroms"],
            allow_missing=True,
        ),
        pops="relate/{panel}-{pops}.poplabels",
    output:
        expand(
            [
                "relate/{panel}-{pops}-popsize_chr{chr}.anc.gz",
                "relate/{panel}-{pops}-popsize_chr{chr}.mut.gz",
                "relate/{panel}-{pops}-popsize_chr{chr}.dist.gz",
            ],
            chr=config["chroms"],
            allow_missing=True,
        ),
        "relate/{panel}-{pops}-popsize_avg.rate",
        "relate/{panel}-{pops}-popsize.chr",
        "relate/{panel}-{pops}-popsize.coal",
        "relate/{panel}-{pops}-popsize.pairwise.bin",
        "relate/{panel}-{pops}-popsize.pairwise.coal",
        "relate/{panel}-{pops}-popsize.pdf",
    log:
        "relate/{panel}-{pops}-popsize.log",
    params:
        input="relate/{panel}-{pops}",
        output="relate/{panel}-{pops}-popsize",
        first_chr=min(config["chroms"]),
        last_chr=max(config["chroms"]),
        dist="relate/{panel}-{pops}-popsize_chr*.dist",
    threads: MAX_THREADS
    shell:
        "{config[relate][path]}/scripts/EstimatePopulationSize/EstimatePopulationSize.sh"
        " --input {params.input}"
        " --output {params.output}"
        " --mu {config[relate][mu]}"
        " --poplabels {input.pops}"
        " --years_per_gen {config[gen_time]}"
        " --first_chr {params.first_chr}"
        " --last_chr {params.last_chr}"
        " --threads {threads}"
        " --seed 122 &> {log} && "
        "parallel -P {threads} 'gzip -c {{}} > {{}}.gz' ::: {params.dist}" # tidy up non-gzipped output


rule relate_map_mutations:
    input:
        anc="relate/{panel}-{pops}-popsize_chr{chr}.anc.gz",
        mut="relate/{panel}-{pops}-popsize_chr{chr}.mut.gz",
        haps="haps/{panel}-{pops}_chr{chr}.haps",
        smpl="haps/{panel}-{pops}_chr{chr}.sample",
    output:
        anc="relate/{panel}-{pops}-popsize-allsnps_chr{chr}.anc.gz",
        mut="relate/{panel}-{pops}-popsize-allsnps_chr{chr}.mut.gz",
        dst="relate/{panel}-{pops}-popsize-allsnps_chr{chr}.dist.gz",
    log:
        "relate/{panel}-{pops}-popsize-allsnps_chr{chr}.log",
    params:
        input=lambda wildcards, input: trim_ext(input[0], 2),
        stub=lambda wildcards, input, output: trim_ext(output[0], 2),
    shell:
        "{config[relate][path]}/bin/RelateExtract"
        " --mode MapMutations"
        " --anc {input.anc}"
        " --mut {input.mut}"
        " --haps {input.haps}"
        " --sample {input.smpl}"
        " --output {params.stub} &> {log} && "
        "gzip {params.stub}.mut && "
        "gzip {params.stub}.dist && "
        "cp {input.anc} {output.anc}"


rule relate_is_not_mapping_chr:
    input:
        mut="relate/{panel}-{pops}-popsize-allsnps_chr{chr}.mut.gz",
    output:
        temp("relate/{panel}-{pops}-popsize-allsnps_chr{chr}_unmapped.tsv"),
    shell:
        r"""gunzip -c {input.mut} | awk -F ';' '$7==1 {{print "{wildcards.chr}:"$2"\t"$4}}' > {output}"""


checkpoint relate_is_not_mapping:
    input:
        expand(
            "relate/{panel}-{pops}-popsize-allsnps_chr{chr}_unmapped.tsv", chr=config["chroms"], allow_missing=True,
        ),
    output:
        "relate/{panel}-{pops}-popsize-allsnps_unmapped.tsv.gz",
    shell:
        r"printf 'pos\trsid\n'| gzip -c > {output} && cat {input} | gzip -c >> {output}"


rule relate_is_flipped_chr:
    input:
        mut="relate/{panel}-{pops}-popsize-allsnps_chr{chr}.mut.gz",
    output:
        temp("relate/{panel}-{pops}-popsize-allsnps_chr{chr}_flipped.tsv"),
    shell:
        r"""gunzip -c {input.mut} | awk -F ';' '$8==1 {{print "{wildcards.chr}:"$2"\t"$4}}' > {output}"""


checkpoint relate_is_flipped:
    input:
        expand(
            "relate/{panel}-{pops}-popsize-allsnps_chr{chr}_flipped.tsv", chr=config["chroms"], allow_missing=True,
        ),
    output:
        "relate/{panel}-{pops}-popsize-allsnps_flipped.tsv.gz",
    shell:
        r"printf 'pos\trsid\n'| gzip -c > {output} && cat {input} | gzip -c >> {output}"


def relate_sample_branch_lengths_inputs(wildcards):
    """Get the reference assembly and chromosome for the given wildcards"""
    reference = "GRCh37"

    # noinspection PyUnresolvedReferences
    meta_file = checkpoints.variant_metadata.get(reference=reference, rsid=wildcards.rsid).output[0]
    with open(meta_file) as fin:
        meta = json.load(fin)

    params = dict(**wildcards)
    params["chr"] = meta["chrom"]
    params["reference"] = reference

    return {
        "anc": expand("relate/{panel}-{pops}-popsize-allsnps_chr{chr}.anc.gz", **params),
        "mut": expand("relate/{panel}-{pops}-popsize-allsnps_chr{chr}.mut.gz", **params),
        "dist": expand("relate/{panel}-{pops}-popsize-allsnps_chr{chr}.dist.gz", **params),
        "coal": expand("relate/{panel}-{pops}-popsize.coal", **params),
        "meta": expand("variants/metadata/{reference}/{rsid}.json", **params),
    }


rule relate_sample_branch_lengths:
    input:
        unpack(relate_sample_branch_lengths_inputs),
    output:
        anc=temp("relate/{panel}/{pops}/{rsid}/{panel}-{pops}-popsize-allsnps-{rsid}.anc"),
        mut=temp("relate/{panel}/{pops}/{rsid}/{panel}-{pops}-popsize-allsnps-{rsid}.mut"),
        dist=temp("relate/{panel}/{pops}/{rsid}/{panel}-{pops}-popsize-allsnps-{rsid}.dist"),
        timeb="relate/{panel}/{pops}/{rsid}/{panel}-{pops}-popsize-allsnps-{rsid}.timeb",
    log:
        "relate/{panel}/{pops}/{rsid}/{panel}-{pops}-popsize-allsnps-{rsid}.timeb.log",
    params:
        input=lambda wildcards, input: trim_ext(input[0], 2),
        output=lambda wildcards, input, output: trim_ext(output[0]),
        pos=lambda wildcards, input: json.load(open(str(input.meta)))["start"],
    shell:
        "{config[relate][path]}/scripts/SampleBranchLengths/SampleBranchLengths.sh"
        " --input {params.input}"
        " --output {params.output}"
        " --mu {config[relate][mu]}"
        " --coal {input.coal}"
        " --num_samples {config[relate][num_samples]}"
        " --first_bp {params.pos}"
        " --last_bp {params.pos}"
        " --dist {input.dist}"
        " --format {RELATE_FORMAT_BINARY}"
        " --seed 122 &> {log}"
