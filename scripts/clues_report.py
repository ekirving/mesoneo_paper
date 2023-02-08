#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json
import sys

import click
import pandas as pd
from scipy.stats.distributions import chi2

ANCESTRIES = ["ALL", "ANA", "CHG", "WHG", "EHG"]


@click.command()
@click.option("--data", "data_tsv", metavar="<file>", help="SNP data", type=click.Path(exists=True), required=True)
@click.option("--columns", metavar="<col,col>", help="rsID columns", required=True)
@click.option("--info", "info_tsv", metavar="<file>", help="INFO scores", type=click.Path(exists=True), required=True)
@click.option("--dataset", metavar="<string>", help="Name of the dataset", required=True)
@click.option("--population", metavar="<string>", help="Name of the population", required=True)
@click.option("--mode", metavar="<string>", help="Clues mode", required=True)
@click.option("--ancestry", metavar="<string>", help="Ancestral path", type=click.Choice(ANCESTRIES), required=True)
@click.option("--output", metavar="<file>", type=click.Path(writable=True), help="Output filename", required=True)
def clues_report(data_tsv, columns, info_tsv, dataset, population, mode, ancestry, output):
    """
    Generate a CLUES report
    """
    # get the list of SNPs to load
    data = pd.read_table(data_tsv)
    info = pd.read_table(info_tsv)

    snps = set()
    for col in columns.split(","):
        snps.update(set(data[col]))

    print(f"INFO: Loading {len(snps):,} SNPs", file=sys.stderr)

    rows = []

    for rsid in snps:
        # load the SNP metadata
        meta_file = f"variants/metadata/GRCh37/{rsid}.json"
        try:
            with open(meta_file) as fin:
                meta = json.load(fin)

                # extract the phenotypes
                meta["pubmed"] = "; ".join(sorted(set(assoc.get("pubmedid") for assoc in meta.get("gwascat", []))))
                meta["gwascat"] = "; ".join(sorted(set(assoc.get("phenotype") for assoc in meta.get("gwascat", []))))

        except FileNotFoundError:
            print(f"ERROR: Missing metadata file {meta_file}", file=sys.stderr)
            meta = None

        # load the model parameters
        mod_file = f"clues/{dataset}/{population}/{rsid}/{dataset}-{population}-{rsid}-{mode}-{ancestry}-any.json"
        try:
            with open(mod_file) as fin:
                model = json.load(fin)

        except FileNotFoundError:
            print(f"WARNING: Missing model {mod_file}", file=sys.stderr)
            model = None

        if meta and model:
            model.pop("sex", None)

            # extract the first (and only) epoch
            epoch, s = list(model.pop("epochs").items())[0]

            model["epoch"] = epoch
            model["s"] = s

            rows.append({**meta, **model})

    print(f"INFO: Loaded {len(rows):,} models for the report", file=sys.stderr)

    # convert to a df
    df = pd.DataFrame(rows)
    df["chrom"] = df.chrom.astype("int64")
    # convert the log-likelihood ratio into a p-value
    # https://en.wikipedia.org/wiki/Wilks%27_theorem
    df["p.value"] = df["logLR"].apply(lambda logLR: chi2.sf(2 * logLR, 1))

    # avoid duplicate column
    info.drop("rsid", axis=1)

    # merge the INFO scores
    merge = pd.merge(df, info, how="left", on=["chrom", "start"])
    merge = merge.sort_values(by=["chrom", "start"])
    merge.to_csv(output, sep="\t", index=False)

    print("INFO: Finished!", file=sys.stderr)


if __name__ == "__main__":
    clues_report()
