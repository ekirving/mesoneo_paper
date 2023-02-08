#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from snakemake.io import temp

# define the region or SNP positions of each loci
LOCI_SNPs = {
    "lct": ["rs1438307", "rs4988235"],
    "plague": ["rs2549794", "rs1052025", "rs11571319", "rs17473484", "rs2248374"],
}
LOCI_POS = {
    "lct": ["2:136608646", "2:136499166"],
    "plague": ["5:96244549", "18:77287776", "2:204738938", "5:114915460", "5:96235896"],
}


rule extract_genotypes_shotgun:
    """
    Extract the genotypes for a locus
    """
    input:
        vcf=lambda wildcards: config["samples"][wildcards.dataset]["genotypes"],
    output:
        tsv="binned/{dataset}-{locus}-genotypes.tsv",
    params:
        locus=lambda wildcards: LOCI_POS[wildcards.locus],
    shell:
        r"bcftools view {input.vcf} {params.locus} | "
        r"bcftools query --format '[%ID\t%REF\t%ALT\t%INFO/AA\t%SAMPLE\t%GT\n]\n' | "
        r"sed 's/|||//' | tr [acgt] [ACGT] > {output.tsv}"


rule plot_binned_freqs_shotgun:
    """
    Plot the binned frequencies for a locus from the shotgun dataset
    """
    input:
        tsv=config["samples"]["ancestral_paths_new"]["metadata"],
        imp="binned/ancestral_paths_new-{locus}-genotypes.tsv",
        lik="binned/neo_likelihoods-{locus}-genotypes.tsv",
    output:
        png="binned/{locus}-binned-calls-shotgun.png",
    shell:
        "Rscript scripts/binned_frequencies_shotgun.R"
        " --samples {input.tsv}"
        " --imputed {input.imp}"
        " --likelihoods {input.lik}"
        " --output {output.png}"


rule download_1240k_data:
    """
    Download the public 1240k array data

    https://reich.hms.harvard.edu/allen-ancient-dna-resource-aadr-downloadable-genotypes-present-day-and-ancient-dna-data
    """
    output:
        anno="data/1240k/v52.2_1240K_public.anno",
        geno="data/1240k/v52.2_1240K_public.geno",
        ind="data/1240k/v52.2_1240K_public.ind",
        snp="data/1240k/v52.2_1240K_public.snp",
        tar=temp("data/1240k/v52.2_1240K_public.tar"),
    shell:
        "wget --quiet -P data/1240k/ https://reichdata.hms.harvard.edu/pub/datasets/amh_repo/curated_releases/V52/V52.2/SHARE/public.dir/v52.2_1240K_public.tar && "
        "tar -xf {output.tar} -C data/1240k/"


rule convert_1240k_plink:
    """
    Convert the 1240k data from Eigenstrat format to PLINK format
    """
    input:
        geno="data/1240k/v52.2_1240K_public.geno",
        ind="data/1240k/v52.2_1240K_public.ind",
        snp="data/1240k/v52.2_1240K_public.snp",
    output:
        bed="data/1240k/plink/v52.2_1240K_public.bed",
        bim="data/1240k/plink/v52.2_1240K_public.bim",
        fam="data/1240k/plink/v52.2_1240K_public.fam",
        par="data/1240k/plink/convert.par",
    log:
        log="data/1240k/plink/convert.log",
    shell:
        r"printf '"
        r"genotypename:    {input.geno}\n"
        r"snpname:         {input.snp}\n"
        r"indivname:       {input.ind}\n"
        r"outputformat:    PACKEDPED\n"
        r"genotypeoutname: {output.bed}\n"
        r"snpoutname:      {output.bim}\n"
        r"indivoutname:    {output.fam}\n"
        r"familynames:     NO\n' > {output.par} && "
        r"convertf -p {output.par} &> {log}"


rule extract_genotypes_1240k:
    """
    Extract the genotypes for a locus from the 1240k dataset
    """
    input:
        bed="data/1240k/plink/v52.2_1240K_public.bed",
        bim="data/1240k/plink/v52.2_1240K_public.bim",
        fam="data/1240k/plink/v52.2_1240K_public.fam",
    output:
        snps=temp("binned/1240k-{locus}.snps"),
        ped=temp("binned/1240k-{locus}.ped"),
        map=temp("binned/1240k-{locus}.map"),
        sex=temp("binned/1240k-{locus}.nosex"),
        tsv="binned/1240k-{locus}.tsv",
    log:
        log="binned/1240k-{locus}.log",
    params:
        snps=lambda wildcards: r"\n".join(LOCI_SNPs[wildcards.locus]),
        prefix="data/1240k/plink/v52.2_1240K_public",
        output="binned/1240k-{locus}",
        colnames=lambda wildcards: r"\t".join(
            ["fid", "sampleId", "father", "mother", "sex", "phenotype"] + LOCI_SNPs[wildcards.locus]
        ),
    shell:
        r"printf '{params.snps}' > {output.snps} && "
        r"plink --bfile {params.prefix} --extract {output.snps} --recode tab --out {params.output} &> {log} && "
        r"printf '{params.colnames}\n' > {output.tsv} && "
        r"cat {output.ped} >> {output.tsv}"


rule plot_binned_freqs_1240k:
    """
    Plot the binned frequencies for a locus from the 1240k dataset
    """
    input:
        tsv="data/1240k/Le_et_al_2022_Table_S1.txt",
        geno="binned/1240k-{locus}.tsv",
    output:
        png="binned/{locus}-binned-calls-1240k.png",
    shell:
        "Rscript scripts/binned_frequencies_1240k.R"
        " --samples {input.tsv}"
        " --geno {input.geno}"
        " --output {output.png}"
