#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

from urllib.parse import quote_plus

import click
import pandas as pd
import requests

EBI_OLS_API = "https://www.ebi.ac.uk/ols/api/ontologies/efo/terms/{}/ancestors"


@click.command()
@click.option("--tsv", "tsv_file", metavar="<file>", help="TSV file", type=click.Path(exists=True), required=True)
@click.option("--output", metavar="<file>", type=click.Path(writable=True), help="Output filename", required=True)
def gwascat_ontology(tsv_file, output):
    """
    Annotate the GWAS Catalog with the full paths of the trait ontology
    """
    gwascat = pd.read_csv(tsv_file, sep="\t", low_memory=False)

    # only use the first of any multiple trait URLs (multiple URLs are used when an older trait has become obsolete)
    gwascat["MAPPED_TRAIT_URI"] = gwascat["MAPPED_TRAIT_URI"].str.extract(pat=r"^(http:[\w.\/]+)")

    # get the unique list of traits
    traits = gwascat["MAPPED_TRAIT_URI"].dropna().unique()

    ontology = []

    for trait_url in traits:
        # the EBI API expects the URL to be double encoded
        uri = quote_plus(quote_plus(trait_url))
        r = requests.get(EBI_OLS_API.format(uri), headers={"Content-Type": "application/json"})

        if not r.ok:
            r.raise_for_status()

        data = r.json()

        labels = []

        for term in data["_embedded"]["terms"]:
            if term["obo_id"] and term["obo_id"].startswith("EFO"):
                labels.append(term["label"])

        ontology.append({"MAPPED_TRAIT_URI": trait_url, "ONTOLOGY_LABELS": ";".join(labels)})

    # convert to df
    ontology = pd.DataFrame(ontology)

    # left join onto the original df
    gwascat = pd.merge(gwascat, ontology, on="MAPPED_TRAIT_URI", how="left")

    cols = {
        "SNPS": "rsid",
        "STRONGEST SNP-RISK ALLELE": "GwasCat_RA",
        "PUBMEDID": "GwasCat_PubMedID",
        "MAPPED_GENE": "GwasCat_Genes",
        "RISK ALLELE FREQUENCY": "GwasCat_RAF",
        "P-VALUE": "GwasCat_P",
        "OR or BETA": "GwasCat_OrorBeta",
        "DISEASE/TRAIT": "GwasCat_Phenotype",
        "INITIAL SAMPLE SIZE": "GwasCat_Samples",
        "ONTOLOGY_LABELS": "GwasCat_Ontology",
    }

    gwascat = gwascat[cols.keys()]
    gwascat.rename(columns=cols, inplace=True)

    # only use the first of any multiple rsIDs
    gwascat["rsid"] = gwascat["rsid"].str.extract(pat=r"^(rs\d+)")
    gwascat["GwasCat_RA"] = gwascat["GwasCat_RA"].str.extract(pat=r"^rs\d+-(.)")

    # drop any GWAS hits where the RA is not listed
    gwascat = gwascat[gwascat["GwasCat_RA"].isin(["A", "C", "G", "T", "?"])]

    # save to disk
    gwascat.to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    gwascat_ontology()
