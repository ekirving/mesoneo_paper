#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json
import os
import sys
from collections import defaultdict
from pprint import pprint

import click
import pysam
import yaml

sys.path.append(os.getcwd())

from scripts.utils import get_samples

# ancestral path codes
PATH_ANA = "1"  # Anatolian Farmers         -> Neolithic
PATH_CHG = "2"  # Caucasus Hunter-gatherers -> Yamnaya
PATH_WHG = "3"  # Western Hunter-gatherers  -> Neolithic
PATH_EHG = "4"  # Eastern Hunter-gatherers  -> Yamnaya
PATH_EHG_WHG = "7"  # North European ancestry (WHG or EHG path) but unable to be more specific
PATH_ANA_CHG = "8"  # West Asian ancestry (CHG or Anatolian path) but unable to be more specific
PATH_UNKNOWN = "0"  # Unable to assign specific path (This labels 0,5,6,9 and 10)


@click.command()
@click.option("--vcf", "vcf_file", metavar="<file>", help="VCF file", type=click.Path(exists=True), required=True)
@click.option("--chr", "chrom", metavar="<chr>", help="Chromosome of the variant", required=True)
@click.option("--pos", metavar="<int>", help="Position of the variant", type=int, required=True)
@click.option("--dataset", metavar="<string>", help="Name of the dataset", required=True)
@click.option("--population", metavar="<string>", help="Name of the population", required=True)
@click.option("--output", metavar="<file>", type=click.File("w"), help="Output filename", required=True)
def ancestral_frequencies(vcf_file, chrom, pos, dataset, population, output):
    """
    Extract the SNP frequencies for each of the ancestral paths
    """
    with open("config.yaml") as fin:
        config = yaml.safe_load(fin)

    # get all the samples in the current analysis group
    samples = get_samples(config, dataset, population)

    # load the VCF file with the sample genotypes
    vcf = pysam.VariantFile(vcf_file)

    try:
        # fetch the record from the VCF
        rec = next(vcf.fetch(chrom, pos - 1, pos))

    except (StopIteration, ValueError):
        # variant not in the VCF
        raise RuntimeError("SNP {}:{} not found in {}".format(chrom, pos, vcf_file))

    data = {
        "ancient": defaultdict(dict),
        "modern": defaultdict(dict),
    }

    for sample, sample_row in samples.iterrows():
        key = "modern" if sample_row["age"] == 0 else "ancient"

        try:
            rec.samples[sample]
        except KeyError:
            print("WARNING: Sample {} not found in the VCF".format(sample), file=sys.stderr)
            continue

        for geno, path in zip(rec.samples[sample].alleles, rec.samples[sample].get("AP", "")[0].split("|")):
            # count the ancestry specific genotypes
            data[key][path][geno] = data[key][path].get(geno, 0) + 1

    pprint(data)

    json.dump(data, output, indent=2, sort_keys=True)


if __name__ == "__main__":
    ancestral_frequencies()
