#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import click
import pandas as pd


@click.command()
@click.option(
    "--report", "report_file", metavar="<file>", help="Phenotypes", type=click.Path(exists=True), required=True
)
@click.option("--ancestral", "anc_file", metavar="<file>", help="AA file", type=click.Path(exists=True), required=True)
@click.option("--output", metavar="<file>", type=click.File("w"), help="Output filename", required=True)
def fix_derived(report_file, anc_file, output):
    """
    Fix the bad derived allele calls
    """
    report_df = pd.read_table(report_file)
    anc = pd.read_table(anc_file)

    # e.g. 1:13011:T:G
    anc["chrom"], anc["start"], anc["ref"], anc["alt"] = anc["variant"].str.split(":").str

    # cast to float64
    anc["start"] = anc.start.astype("float64")

    # join the REF/ALT calls
    report_df = pd.merge(report_df, anc, on=["chrom", "start"])

    # get the bad values
    report_df = report_df[(report_df["derived"] != report_df["ref"]) & (report_df["derived"] != report_df["alt"])]

    report_df.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    fix_derived()
