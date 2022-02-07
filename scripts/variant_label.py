#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json
import re
import textwrap
from collections import defaultdict, OrderedDict

import click
import pysam

LABEL_PMID_LENGTH = 10  # 8 for the PMID + 2 for the delimiter
LABEL_MAX_LENGTH = 58


@click.command()
@click.option("--vcf", "vcf_file", metavar="<file>", type=click.Path(exists=True), help="VCF file", required=True)
@click.option("--meta", "meta_file", metavar="<file>", type=click.File("r"), help="Metadata JSON file", required=True)
@click.option("--out", "output_file", metavar="<file>", type=click.File("w"), help="Output JSON file", required=True)
def variant_label(vcf_file, meta_file, output_file):
    """
    Parse rsID metadata to make a label for the selection plots
    """
    metadata = json.load(meta_file)

    # load the VCF file with the sample genotypes
    vcf = pysam.VariantFile(vcf_file)

    try:
        # fetch the record from the VCF
        rec = next(vcf.fetch(metadata["chrom"], int(metadata["start"]) - 1, int(metadata["end"])))

    except (StopIteration, ValueError):
        # variant not in the VCF
        raise RuntimeError("SNP {}:{} not found in {}".format(metadata["chrom"], metadata["start"], vcf_file))

    if metadata.get("genes", "") == "":
        metadata["genes"] = "N/A"
    else:
        # deduplicate, but keep the original order
        genes = re.split("[,;]+", metadata["genes"])
        metadata["genes"] = "; ".join(list(OrderedDict.fromkeys(genes)))

    # compose the title for the rsID
    title = "{rsid} | chr{chrom}:{start} | Gene(s): {genes} | {ancestral}/{derived}".format(**metadata)

    group = defaultdict(dict)

    # group the GWAS associations by risk allele and phenotype
    for record in metadata.get("gwascat", []):
        if record["phenotype"] not in group[record["allele"]]:
            group[record["allele"]][record["phenotype"]] = [record["pubmedid"]]
        else:
            group[record["allele"]][record["phenotype"]].append(record["pubmedid"])

    gwascat = []

    # summarize the GWAS Catalog associations
    for allele in group:
        # skip GWAS associations for alleles that are not present in this dataset
        for pheno in group[allele]:
            pmid = "; ".join(sorted(set(group[allele][pheno]), reverse=True))

            # handle edge case of too many PMIDs which obscure the phenotype
            if len(pmid) > LABEL_PMID_LENGTH * 4:
                pmid = textwrap.shorten(pmid, width=LABEL_PMID_LENGTH * 4)

            # handle phenotype labels that are too long
            while len(pheno) + len(pmid) > LABEL_MAX_LENGTH:
                if "(" in pheno:
                    # remove any trailing bracketed terms (as these are often redundant)
                    pheno = "(".join(pheno.split("(")[:-1]).strip()
                else:
                    # truncate the pheno
                    pheno = textwrap.shorten(pheno, width=LABEL_MAX_LENGTH - len(pmid))

            gwascat.append("{allele}: {pheno} (PMID: {pmid})".format(allele=allele, pheno=pheno, pmid=pmid))

    label = {"title": title, "gwascat": sorted(gwascat)}

    # save the label
    json.dump(label, output_file, indent=2)


if __name__ == "__main__":
    variant_label()
