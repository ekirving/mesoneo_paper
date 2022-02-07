#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json
import os
import subprocess
import sys

import click
import pandas as pd
import yaml

sys.path.append(os.getcwd())

from scripts.utils import gen_dict_extract


@click.command()
@click.option("--rsid", metavar="<string>", help="RefSeq ID", required=True)
@click.option("--var", "var_file", metavar="<file>", type=click.File("r"), help="Ensembl VAR record", required=True)
@click.option("--vep", "vep_file", metavar="<file>", type=click.File("r"), help="Ensembl VEP record", required=True)
@click.option("--gwas", "gwas_file", metavar="<file>", type=click.File("r"), help="GWAS Catalog file", required=True)
@click.option("--out", "output_file", metavar="<file>", type=click.File("w"), help="Output file", required=True)
def variant_metadata(rsid, var_file, vep_file, gwas_file, output_file):
    """
    Parse the Ensembl record to retrieve the rsID metadata
    """
    with open("config.yaml") as fin:
        config = yaml.safe_load(fin)

    try:
        var = json.load(var_file)
        vep = json.load(vep_file)
    except json.decoder.JSONDecodeError as e:
        print("ERROR loading {} and {}".format(var_file, vep_file))
        raise e

    # get the clinical significance data from the VEP record
    clin_sig = set()
    clin_sig_allele = set()
    severity = set()
    gwas_records = []

    for val in gen_dict_extract("clin_sig_allele", vep):
        for assoc in val.split(";"):
            allele, significance = assoc.split(":")
            clin_sig.add(significance)
            clin_sig_allele.add(allele)

    for val in gen_dict_extract("most_severe_consequence", vep):
        severity.add(val)

    try:
        # get all the GWAS Catalog records for the target rsID
        gwas_target = pd.read_table(gwas_file, dtype=str).set_index("rsid", drop=False).fillna("").loc[[rsid]]
    except KeyError:
        # no matching rows (i.e. this is a neutral SNP)
        gwas_ra = None
        gwas_genes = ""

    else:
        # handle multiple GWAS associations for this rsID
        for idx, gwascat in gwas_target.iterrows():
            # extract the GWAS Catalog allele
            gwascat_allele = gwascat["GwasCat_RA"].split("-").pop()

            # drop any GWAS studies where the RA is not listed
            if gwascat_allele not in ["A", "C", "G", "T", "?"]:
                continue

            gwas_records.append(
                {
                    "pubmedid": gwascat["GwasCat_PubMedID"],
                    "genes": gwascat["GwasCat_Genes"],
                    "allele": gwascat_allele,
                    "raf": gwascat["GwasCat_RAF"],
                    "p": gwascat["GwasCat_P"],
                    "ororbeta": gwascat["GwasCat_OrorBeta"],
                    "phenotype": gwascat["GwasCat_Phenotype"],
                    "samples": gwascat["GwasCat_Samples"],
                    "ontology": gwascat["GwasCat_Ontology"],
                }
            )

        # get the most common risk allele (some rsIDs have multiple)
        gwas_ra = gwas_target["GwasCat_RA"].value_counts().index[0].split("-").pop().replace("?", "")

        # get all the associated genes
        gwas_genes = "/".join(set(gwas_target["GwasCat_Genes"].tolist()))

    # always use the first mapping
    mappings = var.get("mappings", [])
    mapping = mappings[0] if len(mappings) > 0 else {}

    # get the ancestral and derived alleles
    ancestral = mapping.get("ancestral_allele")

    # load the 1000G VCF and check the alleles
    cmd = (
        r"bcftools view --types snps --max-alleles 2 "
        r" {path}/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5.20130502.genotypes.vcf.gz {chr}:{pos} | "
        r"bcftools query --format '%REF %ALT'".format(
            path=config["1000G"]["vcf_path"], chr=mapping.get("seq_region_name"), pos=mapping.get("start")
        )
    )
    alleles = subprocess.check_output([cmd], shell=True).decode("utf-8").split()

    if alleles == "":
        # use the Ensembl annotation
        alleles = mapping.get("allele_string").split("/")

    # the derived is whichever allele is left
    derived = (set(alleles) - {ancestral}).pop()

    # save the metadata for the variant
    metadata = {
        "rsid": var.get("name"),
        "type": var.get("var_class", var.get("error")),
        "chrom": mapping.get("seq_region_name"),
        "start": mapping.get("start"),
        "end": mapping.get("end"),
        "allele": gwas_ra,
        "genes": gwas_genes,
        "alleles": mapping.get("allele_string"),
        "ancestral": ancestral,
        "derived": derived,
        "minor": var.get("minor_allele"),
        "clin_allele": "/".join(set(clin_sig_allele)),
        "clin_sig": "/".join(set(clin_sig)),
        "severity": "/".join(set(severity)),
        "gwascat": gwas_records,
    }

    # save the metadata
    json.dump(metadata, output_file, indent=2)


if __name__ == "__main__":
    variant_metadata()
