#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Customised version of https://github.com/35ajstern/clues/blob/master/plot_traj.py
"""

__author__ = "Aaron Stern"
__email__ = "ajstern@berkeley.edu"

import argparse
import json
import warnings
from math import ceil

import numpy as np
from matplotlib import pyplot as plt
from scipy.stats.distributions import chi2

parser = argparse.ArgumentParser()
parser.add_argument("inputPrefix", type=str)
parser.add_argument("figurePrefix", type=str)
parser.add_argument("--ext", type=str, default="pdf")
parser.add_argument("--gen-time", type=int, default="28")
parser.add_argument("--params", type=str)
parser.add_argument("--label", type=str)
parser.add_argument("--ancestry", type=str)
parser.add_argument("--sex", type=str)
args = parser.parse_args()

epochs = np.load("{}.epochs.npy".format(args.inputPrefix))
freqs = np.load("{}.freqs.npy".format(args.inputPrefix))
logpost = np.load("{}.post.npy".format(args.inputPrefix))

with open(args.label) as fin:
    label = json.load(fin)

with open(args.params) as fin:
    params = json.load(fin)

f, ax = plt.subplots(1, 1)
f.set_size_inches(20, 10)

xmin = int(min(epochs))
xmax = min(int(max(epochs)), round(13665 / args.gen_time))  # TODO parameterize this

xticks = range(xmin, xmax + 1, round(1000 / args.gen_time))

# flip the x-axis
epochs = epochs * -1
xticks = [tick * -1 for tick in xticks]
xlabels = [-int(tick * args.gen_time / 1000) for tick in xticks]

desc = {
    "ancient": "Ancient samples only",
    "modern": "Modern 1000G data only",
    "both": "Ancient samples plus modern 1000G data",
}

subtitle = desc[params["mode"]]

ancestries = {
    "ALL": "All ancestries",
    "ANA": "Anatolian Farmers",
    "CHG": "Caucasus Hunter-gatherers",
    "WHG": "Western Hunter-gatherers",
    "EHG": "Eastern Hunter-gatherers",
}

subtitle += " | " + ancestries[args.ancestry]

data = []
for epoch, s in params["epochs"].items():
    # convert the log-likelihood ratio into a p-value
    # https://en.wikipedia.org/wiki/Wilks%27_theorem
    params["p.value"] = chi2.sf(2 * params["logLR"], 1)
    data.append(
        "logLR = {:.2f} | p = {:.2e} | epoch = {} | s = {:.5f}".format(params["logLR"], params["p.value"], epoch, s)
    )

subtitle += "\n" + "\n".join(data)

# ignore MatplotlibDeprecation warnings
warnings.filterwarnings("ignore")

plt.pcolormesh(epochs[:-1], freqs, np.exp(logpost)[:, :])
plt.suptitle(label["title"], x=0.6, fontsize=18)
plt.title(subtitle, fontsize=16, pad=10)
plt.axis((-xmax, -xmin, 0, 1.0))
ax.yaxis.tick_right()
ax.yaxis.set_label_position("right")
plt.ylabel("Derived Allele Frequency", fontsize=16, labelpad=20)
plt.xlabel("kyr BP", fontsize=16)
plt.xticks(ticks=xticks, labels=xlabels, fontsize=16)
plt.yticks(fontsize=18)

cbar = plt.colorbar(ax=[ax], location="left")
plt.clim(0, 0.5)
cbar.ax.set_ylabel("Posterior prob.\n\n", rotation=270, fontsize=16, labelpad=-40)
cbar.ax.set_yticklabels([0.0, 0.1, 0.2, 0.3, 0.4, "≥0.5"])
cbar.ax.tick_params(labelsize=18)

# print in two columns
per_column = ceil(len(label["gwascat"]) / 2)

# add GWAS associations to bottom
for i, label in enumerate(label["gwascat"]):
    if i < per_column:
        x, y = 0.15, i * -0.03
        align = "left"
    else:
        x, y = 0.95, (i - per_column) * -0.03
        align = "right"

    plt.figtext(x, y, label, horizontalalignment=align, fontsize=16)

plt.savefig("%s.%s" % (args.figurePrefix, args.ext), format=args.ext, bbox_inches="tight")
