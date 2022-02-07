#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import os
import sys
from math import ceil, floor

import click
import yaml

sys.path.append(os.getcwd())

from scripts.utils import get_samples


@click.command()
@click.option("--dataset", metavar="<string>", help="Name of the dataset", required=True)
@click.option("--population", metavar="<string>", help="Name of the population", required=True)
@click.option("--gen-time", metavar="<int>", help="Years per generation", type=int, required=True)
@click.option("--output", metavar="<file>", type=click.File("w"), help="Output filename", required=True)
def clues_time_bins(dataset, population, gen_time, output):
    """
    Generate a time bins file for `clues` with one epoch, that lasts from the oldest sample to the most recent.

    See https://github.com/35ajstern/clues/
    """
    with open("config.yaml") as fin:
        config = yaml.safe_load(fin)

    # get all the samples in the current analysis group
    samples = get_samples(config, dataset, population)
    epoch_start = floor(min(samples["age"]) / gen_time)
    epoch_end = ceil(max(samples["age"]) / gen_time)

    output.write("{:.1f}\n{:.1f}\n".format(epoch_start, epoch_end))


if __name__ == "__main__":
    clues_time_bins()
