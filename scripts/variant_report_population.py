#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import click
import pandas as pd

VARIANT_MIN_DAF = 0.01


@click.command()
@click.option("--report", metavar="<file>", type=click.File("r"), help="Full variants report", required=True)
@click.option("--freqs", metavar="<file>", type=click.File("r"), help="Allele frequencies report", required=True)
@click.option("--output", metavar="<file>", type=click.File("w"), help="Output file", required=True)
def variant_group_report(report, freqs, output):
    """
    Subset the main report to only include rsIDs that are present and segregating in the current analysis group.
    """
    data = pd.read_table(report, na_values="-", low_memory=False)
    freqs = pd.read_table(freqs, low_memory=False)
    freqs["AF"] = freqs["AC"] / freqs["AN"]

    # drop missing data
    data = data[data["start"].notnull()]

    # normalise the column types so we can merge
    data.chrom = data.chrom.astype(str)
    freqs.chr = freqs.chr.astype(str)

    # join the SNP frequencies
    data = pd.merge(data, freqs, left_on=["chrom", "start"], right_on=["chr", "pos"], how="inner", suffixes=("", "_y"))

    # we only want SNPs above a minimum AF
    data = data[data["type"] == "SNP"]
    data = data[data["AF"] > VARIANT_MIN_DAF]

    # BED is zero-based and half open
    data["stop"] = data["start"]
    data["start"] = data["start"] - 1

    # cast to int
    data.start = data.start.astype(int)
    data.stop = data.stop.astype(int)

    # drop all sites with unknown ancestral allele
    data = data[data["ancestral"].isin(["A", "C", "G", "T"])]

    # fill empty records
    data.fillna("", inplace=True)

    # subset columns and sort
    data = data[["chrom", "start", "stop", "rsid", "allele", "ancestral"]]
    data = data.sort_values(["chrom", "start"])

    data.to_csv(output, sep="\t", index=False, header=False)


if __name__ == "__main__":
    variant_group_report()
