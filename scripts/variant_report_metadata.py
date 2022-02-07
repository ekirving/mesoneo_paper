#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import csv
import json
from collections import OrderedDict

import click


@click.command()
@click.option("--ref", "reference", metavar="<string>", help="Reference assembly", required=True)
@click.option("--rsids", "rsid_file", metavar="<file>", type=click.File("r"), help="Input rsIDs file", required=True)
@click.option("--out", "output_file", metavar="<file>", type=click.File("w"), help="Output file", required=True)
def variant_report_metadata(reference, rsid_file, output_file):
    """
    Get the variant metadata from the JSON files.
    """
    variants = []

    for rsid in rsid_file:
        if rsid.startswith("rs"):
            with open("variants/metadata/{}/{}.json".format(reference, rsid.strip())) as fin:
                var = json.load(fin)
                var.pop("gwascat", None)
                variants.append(var)

    # make these columns first
    first = ["rsid", "type", "chrom", "start", "end", "allele"]
    columns = list(OrderedDict.fromkeys(first + list(variants[0].keys())))

    # save the variants to a CSV
    w = csv.DictWriter(output_file, columns, delimiter="\t")
    w.writeheader()
    w.writerows(variants)


if __name__ == "__main__":
    variant_report_metadata()
