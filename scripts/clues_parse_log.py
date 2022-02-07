#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import json

import click

SEXES = ["XX", "XY", "any"]


@click.command()
@click.option("--rsid", metavar="<string>", help="RefSeq ID", required=True)
@click.option("--mode", metavar="<string>", help="Clues mode", required=True)
@click.option("--ancestry", metavar="<string>", help="Ancestral path", required=True)
@click.option("--sex", metavar="<string>", help="Sample sex", type=click.Choice(SEXES), required=True)
@click.option("--log", "log_file", metavar="<file>", type=click.Path(writable=True), help="Log file", required=True)
@click.option("--out", "output", metavar="<file>", type=click.File("w"), help="Output filename", required=True)
def clues_parse_log(rsid, mode, ancestry, sex, log_file, output):
    """
    Parse the Clues log file to extract the information we want.
    """
    with open(log_file) as fin:
        epochs = dict()

        for line in fin:
            if "logLR" in line:
                lnl_ratio = line.split().pop()
            elif "epoch" in line and "selection" in line:
                # handle multiple epochs
                while True:
                    try:
                        epoch, s = next(fin).split()
                        epochs[epoch] = float(s)
                    except ValueError:
                        break

        data = {
            "rsid": rsid,
            "mode": mode,
            "ancestry": ancestry,
            "sex": sex,
            "logLR": float(lnl_ratio),
            "epochs": epochs,
        }

    json.dump(data, output, indent=2)


if __name__ == "__main__":
    clues_parse_log()
