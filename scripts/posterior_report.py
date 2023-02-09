#!/usr/bin/env python
# -*- coding: utf-8 -*-

__author__ = "Evan K. Irving-Pease"
__copyright__ = "Copyright 2020, University of Copenhagen"
__email__ = "evan.irvingpease@gmail.com"
__license__ = "MIT"

import csv
import itertools
import multiprocessing as mp
import warnings

import click
import numpy as np
import ot
import pandas as pd
from scipy.spatial.distance import jensenshannon


def posterior_distance(dataset1, population1, dataset2, population2, snp_type, rsid):
    """
    Calculate the Earth-mover distance (EMD) and Jensen-Shannon distance between each epoch in the posteriors
    """
    try:
        # load the posteriors for both models
        logpost1 = np.load(
            f"clues/{dataset1}/{population1}/{rsid}/{dataset1}-{population1}-{rsid}-ancient-ALL-any.post.npy"
        )
        logpost2 = np.load(
            f"clues/{dataset2}/{population2}/{rsid}/{dataset2}-{population2}-{rsid}-ancient-ALL-any.post.npy"
        )

    except FileNotFoundError as warning:
        # report missing models as a warning
        print(warning)
        return {}

    # convert back from log-posteriors
    posterior1 = np.exp(logpost1)
    posterior2 = np.exp(logpost2)

    # number of frequency discretisation bins
    num_freqs = posterior1.shape[0]

    # the trajectories may not be the same length, so only compare overlapping epochs
    num_epochs = posterior1.shape[1] if posterior1.shape < posterior2.shape else posterior2.shape[1]

    # make a cost matrix with a linear transport distance penalty
    cost_matrix = [[abs(i - j) for i in range(num_freqs)] for j in range(num_freqs)]

    # calculate the earth-mover distance for each epoch in the model
    emd = [
        ot.emd2(
            # correct any numerical overflow issues (i.e. sum(p)>1)
            posterior1[:, epoch] / np.sum(posterior1[:, epoch]),
            posterior2[:, epoch] / np.sum(posterior2[:, epoch]),
            cost_matrix,
        )
        for epoch in range(num_epochs)
    ]

    # ignore warnings, as we skip NA values with np.nansum()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        # calculate the Jensen-Shannon distance between each epoch in the posteriors
        # https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.jensenshannon.html
        js = [
            jensenshannon(
                # correct any numerical underflow issues (i.e. p=0)
                posterior1[:, epoch] + np.finfo(np.float64).eps,
                posterior2[:, epoch] + np.finfo(np.float64).eps,
            )
            for epoch in range(num_epochs)
        ]

    return {
        "type": snp_type,
        "rsid": rsid,
        "epochs": num_epochs,
        "emd_mean": np.nansum(emd) / num_epochs,
        "emd_tss": np.nansum(np.square(emd)),
        "js_mean": np.nansum(js) / num_epochs,
        "js_tss": np.nansum(np.square(js)),
    }


@click.command()
@click.option("--dataset1", metavar="<string>", help="Name of the first dataset", required=True)
@click.option("--population1", metavar="<string>", help="Name of the first population", required=True)
@click.option("--dataset2", metavar="<string>", help="Name of the second dataset", required=True)
@click.option("--population2", metavar="<string>", help="Name of the second population", required=True)
@click.option("--pairs", "pairs_tsv", metavar="<file>", help="Pairings", type=click.Path(exists=True), required=True)
@click.option("--out", "output_file", metavar="<file>", type=click.File("w"), help="Output file", required=True)
def posterior_report(dataset1, population1, dataset2, population2, pairs_tsv, output_file):
    """
    Compare the allele trajectory posteriors from two datasets.
    """
    pairs = pd.read_table(pairs_tsv)

    params = []
    for snp_type in ["gwas", "neutral"]:
        params += list(
            zip(
                itertools.repeat(dataset1),
                itertools.repeat(population1),
                itertools.repeat(dataset2),
                itertools.repeat(population2),
                itertools.repeat(snp_type),
                pairs[snp_type].tolist(),
            )
        )

    with mp.Pool(processes=mp.cpu_count()) as pool:
        data = pool.starmap(posterior_distance, params)

    # drop blank rows
    data = [row for row in data if row]

    # save the report
    w = csv.DictWriter(output_file, data[0].keys(), delimiter="\t", restval="")
    w.writeheader()
    w.writerows(data)


if __name__ == "__main__":
    posterior_report()
