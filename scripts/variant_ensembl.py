#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json
import time

import click
import requests

# number of seconds to wait before making a new request
ENSEMBL_WAIT_TIME = 1

# how many seconds to wait between retries (weighted by number of attempts)
ENSEMBL_RETRY_WAIT = 2

# maximum number of times to retry fetching an Ensembl record before giving up
ENSEMBL_MAX_RETRY = 3

ENSEMBL_VAR = "ensembl_var"
ENSEMBL_VEP = "ensembl_vep"


def query_ensembl_api(reference, rsid, json_file, mode=ENSEMBL_VAR, attempts=0):
    """
    Query ensembl to get the variant record for this rsid
    """
    if reference == "GRCh37":
        server = "http://grch37.rest.ensembl.org"
    elif reference == "GRCh38":
        server = "https://rest.ensembl.org"
    else:
        raise RuntimeError("Unknown reference {}".format(reference))

    if mode == ENSEMBL_VAR:
        url = "/variation/human/{}".format(rsid)
    elif mode == ENSEMBL_VEP:
        url = "/vep/human/id/{}".format(rsid)
    else:
        raise RuntimeError("Unknown service mode {}".format(mode))

    r = requests.get(server + url, headers={"Content-Type": "application/json"})

    if not r.ok:
        # handle Too Many Requests error
        if attempts < ENSEMBL_MAX_RETRY:
            time.sleep(ENSEMBL_RETRY_WAIT * attempts)
            query_ensembl_api(reference, rsid, json_file, mode, attempts + 1)
            return

        elif r.status_code == 400:
            # some rsIDs might not exist in this particular assembly
            pass
        else:
            r.raise_for_status()

    data = r.json()

    # some rsIDs are aliases for another, so the canonical name may not match our rsID
    if mode == ENSEMBL_VAR and data.get("name") != rsid:
        if data.get("name"):
            data["synonyms"].append(data.get("name"))
            data["synonyms"].remove(rsid)
        data["name"] = rsid

    # save the data
    json.dump(data, json_file, indent=2, sort_keys=True)


@click.command()
@click.option("--ref", metavar="<string>", help="Reference assembly", required=True)
@click.option("--rsid", metavar="<string>", help="RefSeq ID", required=True)
@click.option("--var", "var_file", metavar="<file>", type=click.File("w"), help="Output for VAR record", required=True)
@click.option("--vep", "vep_file", metavar="<file>", type=click.File("w"), help="Output for VEP record", required=True)
def variant_ensembl(ref, rsid, var_file, vep_file):
    # get both Ensembl records
    query_ensembl_api(ref, rsid, var_file, mode=ENSEMBL_VAR)
    query_ensembl_api(ref, rsid, vep_file, mode=ENSEMBL_VEP)
    time.sleep(ENSEMBL_WAIT_TIME)


if __name__ == "__main__":
    variant_ensembl()
