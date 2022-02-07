#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import sys

import click
import pandas as pd

# the rounding precision to use for binning allele frequencies
ALLELE_FREQ_PRECISION = 2

# use a higher precision for pairing the simulated SNPs
SIMULATED_FREQ_PRECISION = 3


@click.command()
@click.option("--gwas", "gwas_tsv", metavar="<file>", help="GWAS SNPs", type=click.Path(exists=True), required=True)
@click.option("--neutrals", metavar="<file>", help="Neutral SNPs", type=click.Path(exists=True), required=True)
@click.option("--ancestral", "anc_tsv", metavar="<file>", help="Ancestral", type=click.Path(exists=True), required=True)
@click.option("--unmapped", "unm_tsv", metavar="<file>", help="Unmapped", type=click.Path(exists=True), required=True)
@click.option("--simulated", is_flag=True)
@click.option("--output", metavar="<file>", type=click.Path(writable=True), help="Output filename", required=True)
@click.option("--ratio", metavar="<int>", help="Ratio of neutral to GWAS", type=int, default=1)
def neutral_pair_snps(gwas_tsv, neutrals, anc_tsv, unm_tsv, simulated, output, ratio):
    """
    Pair GWAS SNPs with neutral SNPs of the same modern frequency.
    """
    precision = SIMULATED_FREQ_PRECISION if simulated else ALLELE_FREQ_PRECISION

    # load the GWAS SNPs
    gwas = pd.read_table(gwas_tsv)
    gwas["AF"] = round(gwas["AC"] / gwas["AN"], precision)

    # load the neutral SNPs
    neut = pd.read_table(neutrals)
    neut["AF"] = round(neut["AC"] / neut["AN"], precision)

    # split SNPs with multiple rsIDs
    gwas["rsid"] = gwas["rsid"].str.split(";")
    gwas = gwas.explode("rsid")

    # drop any SNPs that don't have rsIDs
    gwas = gwas[gwas["rsid"].str.startswith("rs", na=False)]

    # load the ancestral alleles
    ancestral = pd.read_table(anc_tsv)

    # join the ancestral calls, and polarize against the ancestral
    gwas = pd.merge(gwas, ancestral, on="variant")
    gwas.loc[gwas.alt == gwas.ancestral, "AF"] = 1 - gwas.loc[gwas.alt == gwas.ancestral, "AF"]

    # drop gwas SNPs where the ancestral is not in REF/ALT
    gwas = gwas[(gwas.ref == gwas.ancestral) | (gwas.alt == gwas.ancestral)]

    if not simulated:
        # only use the first of any multiple rsIDs
        neut["rsid"] = neut["rsid"].str.extract(pat=r"^(rs\d+)")

        # drop any SNPs that don't have rsIDs
        neut = neut[neut["rsid"].str.startswith("rs", na=False)]

        # join the ancestral calls
        neut = pd.merge(neut, ancestral, on="variant")

        # drop neutral SNPs where the ancestral is not in REF/ALT
        neut = neut[(neut.ref == neut.ancestral) | (neut.alt == neut.ancestral)]

        # polarize against the ancestral
        neut.loc[neut.alt == neut.ancestral, "AF"] = 1 - neut.loc[neut.alt == neut.ancestral, "AF"]

    # apply a random shuffle, so the neutral SNPs are not ordered by position
    neut = neut.sample(frac=1, random_state=122).reset_index(drop=True)

    # don't include neutral SNPs with no modern observations
    neut = neut[neut["AC"] != 0]

    if not simulated:
        # extract just the chr:pos portion of the variant name
        neut["pos"] = neut["variant"].str.extract(pat=r"^(\d+:\d+)")

        # load the list of SNPs that Relate could not map
        unmapped = pd.read_table(unm_tsv)

        # drop all the non-mapping SNPs
        neut = neut[~neut.pos.isin(unmapped["pos"])]

    # bin all the SNPs by DAF
    gwas = gwas.groupby(["chr", "AF"])["rsid"].apply(";".join).reset_index()
    gwas = gwas.rename(columns={"rsid": "gwas"})

    neut = neut.groupby(["chr", "AF"])["rsid"].apply(";".join).reset_index()
    neut = neut.rename(columns={"rsid": "neutrals"})

    # pair the GWAS SNPs with their neutral counterparts
    if not simulated:
        bins = pd.merge(gwas, neut, on=["chr", "AF"], how="left").fillna("")
    else:
        bins = pd.merge(gwas, neut, on=["chr", "AF"], how="inner")

    snps = []
    chrom = ""
    unused = []

    # iterate over the AF bins and draw the same number of SNPs from neutrals
    for idx, af_bin in bins.iterrows():
        # keep track the the unused SNPS (in case we need to borrow one)
        if chrom != af_bin["chr"]:
            chrom = af_bin["chr"]
            unused = []

        pair_gwas = [b for b in af_bin["gwas"].split(";") if b]
        pair_neut = [n for n in af_bin["neutrals"].split(";") if n]

        num_gwas = len(pair_gwas)
        num_neut = len(pair_neut)

        if num_gwas * ratio > num_neut:
            print(
                "WARNING: chr{} - bin {:.2f} - Borrowing {} unpaired neutral SNPs from a lower bin".format(
                    af_bin["chr"], af_bin["AF"], (num_gwas * ratio) - num_neut
                ),
                file=sys.stderr,
            )

        # back fill any missing pairs with unused SNPs from the previous bins
        pair_neut += unused
        unused = pair_neut[(num_gwas * ratio) :]

        # do the pairing
        for gwas_snp in pair_gwas:
            snps.append([af_bin["chr"], gwas_snp, ";".join(pair_neut[:ratio])])
            pair_neut = pair_neut[ratio:]
    # snps.extend(zip(repeat(af_bin["chr"]), pair_gwas, pair_neut))

    # save the pairings
    pairs = pd.DataFrame(snps, columns=["chr", "gwas", "neutral"])
    pairs.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    neutral_pair_snps()
