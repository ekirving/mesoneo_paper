#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import os
import sys
from math import log

import click
import pysam
import yaml

sys.path.append(os.getcwd())

from scripts.utils import get_samples


def get_ancestry_map(dataset):
    """
    Handle the difference in ancestry codes between `ancestral_paths_new` and `ancestral_paths_v3`
    """
    PATH_ANA = "1"  # Anatolian Farmers -> Neolithic
    PATH_CHG = "2"  # Caucasus Hunter-gatherers -> Yamnaya
    PATH_WHG = "3"  # Western Hunter-gatherers -> Neolithic
    PATH_EHG = "4"  # Eastern Hunter-gatherers -> Yamnaya

    ancestry_map = dict()

    if dataset in ["ancestral_paths_new", "chr3_true_paths", "chr3_inferred_paths"]:

        PATH_EHG_WHG = "5"  # North European ancestry (WHG or EHG path) but unable to be more specific
        PATH_ANA_CHG = "6"  # West Asian ancestry (CHG or Anatolian path) but unable to be more specific
        PATH_UNKNOWN = "0"  # Unable to assign specific path (This labels 0,5,6,9 and 10)

        ancestry_map = {
            "ALL": None,
            "ANA": [PATH_ANA, PATH_ANA_CHG],
            "CHG": [PATH_CHG, PATH_ANA_CHG],
            "WHG": [PATH_WHG, PATH_EHG_WHG],
            "EHG": [PATH_EHG, PATH_EHG_WHG],
        }
    elif dataset in ["ancestral_paths_v3", "simulated_relate_painted"]:

        PATH_ANA_BAA = "5"  # Anatolian Farmer -> Bronze Age Anatolian
        PATH_CHG_BAA = "6"  # Caucasus Hunter-gatherers -> Bronze Age Anatolian

        # NB we don't use paths 5 and 6 as they lead to BAA
        ancestry_map = {
            "ALL": None,
            "ANA": [PATH_ANA],
            "CHG": [PATH_CHG],
            "WHG": [PATH_WHG],
            "EHG": [PATH_EHG],
        }

    return ancestry_map


SEXES = ["XX", "XY", "any"]

BASES = ["A", "C", "G", "T"]
BASES_N = BASES + ["N", "0"]

# the minimum number of samples to model the trajectory
MIN_ANCIENT_SAMPLES = 2


@click.command()
@click.option("--vcf", "vcf_file", metavar="<file>", help="VCF file", type=click.Path(exists=True), required=True)
@click.option("--chr", "chrom", metavar="<chr>", help="Chromosome of the variant", required=True)
@click.option("--pos", metavar="<int>", help="Position of the variant", type=int, required=True)
@click.option("--ancestral", metavar="<chr>", help="The ancestral allele", type=click.Choice(BASES_N), required=True)
@click.option("--dataset", metavar="<string>", help="Name of the dataset", required=True)
@click.option("--population", metavar="<string>", help="Name of the population", required=True)
@click.option(
    "--ancestry",
    metavar="<string>",
    help="Ancestry code",
    type=click.Choice(["ALL", "ANA", "CHG", "WHG", "EHG"]),
    required=True,
)
@click.option("--sex", metavar="<string>", help="Sample sex", type=click.Choice(SEXES), required=True)
@click.option("--gen-time", metavar="<int>", help="Years per generation", type=int, required=True)
@click.option("--mod-freq", metavar="<file>", type=click.File("w"), help="Modern frequency filename", required=True)
@click.option("--output", metavar="<file>", type=click.File("w"), help="Output filename", required=True)
def clues_ancient_samples(
    vcf_file, chrom, pos, ancestral, dataset, population, ancestry, sex, gen_time, mod_freq, output
):
    """
    Generate the ancient samples genotype likelihood file for `clues`.

    See https://github.com/35ajstern/clues/
    """
    with open("config.yaml") as fin:
        config = yaml.safe_load(fin)

    # get all the ancient samples in the current analysis group
    samples = get_samples(config, dataset, population)
    ancients = samples[samples["age"] != 0]
    ancients = ancients[ancients["age"].notnull()]
    ancients = ancients.sort_values("age")
    ancestry_map = get_ancestry_map(dataset)

    # also get all the modern samples
    moderns = samples[samples["age"] == 0]

    # apply the sex filter
    if sex != "any":
        ancients = ancients[ancients["sex"] == sex]
        moderns = moderns[moderns["sex"] == sex]

    # load the VCF file with the sample genotypes
    vcf = pysam.VariantFile(vcf_file)

    try:
        # fetch the record from the VCF
        rec = next(vcf.fetch(chrom, pos - 1, pos))

    except (StopIteration, ValueError):
        # variant not in the VCF
        raise RuntimeError(f"SNP {chrom}:{pos} not found in {vcf_file}")

    alleles = [rec.ref] + list(rec.alts)

    if ancestral == "0":
        ancestral = rec.ref

    if len(alleles) > 2:
        raise RuntimeError(f"{chrom}:{pos} {rec.id} SNP is polyallelic {alleles} in {vcf_file}")

    if ancestral == "N":
        raise RuntimeError(f"{chrom}:{pos} {rec.id} Cannot handle SNPs with unknown ancestral allele")

    if ancestral not in alleles:
        raise RuntimeError(f"{chrom}:{pos} {rec.id} Ancestral allele {ancestral} is missing {alleles} in {vcf_file}")

    derived = (set(alleles) - {ancestral}).pop()

    num_samples = 0

    # is the current dataset using genotype likelihoods for ancient samples
    is_likelihood = config["samples"][dataset]["is_likelihood"]

    for sample, sample_row in ancients.iterrows():
        gen = sample_row["age"] / gen_time

        if None in rec.samples[sample].alleles:
            # skip sites without diploid coverage
            continue

        if is_likelihood:
            # get the Phred-scaled genotype likelihoods
            pred_likelihoods = rec.samples[sample].get("PL", None)

            # convert from Phred back into a probability
            gp_diploid = [10 ** (-Q / 10) for Q in pred_likelihoods]

        else:
            # get the diploid genotype probabilities
            gp_diploid = rec.samples[sample].get("GP", None)

            # handle samples without a genotype probability
            if gp_diploid is None:
                # convert regular diploid genotypes into a GP tuple
                alts = rec.samples[sample].alleles.count(rec.alts[0])
                gp_diploid = [0, 0, 0]
                gp_diploid[alts] = 1

        # treat call as diploid when we're not conditioning on ancestry
        if ancestry == "ALL":
            if rec.ref != ancestral:
                # polarise the probabilities
                gp_diploid = reversed(gp_diploid)

            # convert GP calls into pseudo-likelihoods
            geno_ll = [log(geno, 10) if geno > 0 else float("-inf") for geno in gp_diploid]

            # output the pseudo-likelihoods
            output.write("{:f} {:f} {:f} {:f}\n".format(gen, *geno_ll))

            num_samples += 1

        else:
            # get the genotype tuple
            gt = rec.samples[sample].get("GT")

            # treat each call as pseudo-haploid
            for geno, path in zip(gt, rec.samples[sample].get("AP", "")[0].split("|")):
                if path in ancestry_map[ancestry]:
                    if gt == (1, 1):
                        gp_haploid = (gp_diploid[0] + (gp_diploid[1] / 2), gp_diploid[2])
                    elif gt == (0, 0):
                        gp_haploid = (gp_diploid[0], (gp_diploid[1] / 2) + gp_diploid[2])
                    else:
                        # the haploid probability of a het call depend on which call it is
                        if geno == 1:
                            gp_haploid = (gp_diploid[0], gp_diploid[1] + gp_diploid[2])
                        else:
                            gp_haploid = (gp_diploid[0] + gp_diploid[1], gp_diploid[2])

                    if rec.ref != ancestral:
                        # polarise the probabilities
                        gp_haploid = reversed(gp_haploid)

                    # convert GP calls into pseudo-likelihoods
                    geno_ll = [log(geno, 10) if geno > 0 else float("-inf") for geno in gp_haploid]

                    # output the pseudo-likelihoods
                    output.write("{:f} {:f} {:f}\n".format(gen, *geno_ll))

                    num_samples += 1

    if num_samples < MIN_ANCIENT_SAMPLES:
        # output some null records so CLUES doesn't throw an error
        for _ in range(MIN_ANCIENT_SAMPLES - num_samples):
            if ancestry == "ALL":
                output.write("1.0 0.0 -inf -inf\n")
            else:
                output.write("1.0 0.0 -inf\n")

    # calculate the modern frequency
    focal, total = 0, 0

    for sample, sample_row in moderns.iterrows():
        if None in rec.samples[sample].alleles:
            # skip sites without diploid coverage
            continue

        if ancestry == "ALL":
            # count diploid occurrences of the derived allele
            focal += rec.samples[sample].alleles.count(derived)
            total += len(rec.samples[sample].alleles)

        else:
            # count haploid occurrences of the derived allele
            calls = []

            # filter for genotypes belonging to this ancestry
            for call, path in zip(rec.samples[sample].alleles, rec.samples[sample].get("AP", "")[0].split("|")):
                if path in ancestry_map[ancestry]:
                    # count the ancestry specific genotypes
                    calls.append(call)

            focal += calls.count(derived)
            total += len(calls)

    mod_freq.write("{:.4f}".format((focal / total) if total else 0))


if __name__ == "__main__":
    clues_ancient_samples()
