#!/usr/bin/env python3
"""
1.07_appraise_percent_binned_plot.py

Percentage of SingleM mOTU coverage in each sample that is assigned to a MAG
or a GTDB R10-RS226 genome (singlem appraise, 1.06), at species and genus
level, for sediment and seawater, all reads and 10 M read pairs.

One point per sample; black bars mark the pooled value over all samples (the
"total" row written by singlem appraise). The "total" and "average" rows and
"#" comment lines are dropped from the per-sample values.

Input files are named <NN>_<hab>_<reads>_<ref>.<rank>_appraise.tsv, e.g.
00_sed_allreads_mags.species_appraise.tsv. Habitat is taken from the first
digit of NN (0 = sediment, 1 = seawater); a warning is printed if the habitat
in the file name disagrees.

Usage : python3 scripts/1.07_appraise_percent_binned_plot.py <results_dir> <out_prefix>
Output: <out_prefix>.png, <out_prefix>.svg, <out_prefix>.tsv (per-sample values)
"""
import glob
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

results_dir, out_prefix = sys.argv[1], sys.argv[2]

HABITAT  = {"0": "Sediment", "1": "Seawater"}
HAB_NAME = {"sed": "Sediment", "sea": "Seawater"}
READS    = {"allreads": "All reads", "10Mreads": "10 M read pairs"}
REF      = {"mags": "MAGs", "tkrs226genomes": "GTDB R10-RS226"}
COLOURS  = {"Sediment": "#33a02c", "Seawater": "#1f78b4"}
name_re  = re.compile(r"(\d)(\d)_(sed|sea)_(allreads|10Mreads)_(mags|tkrs226genomes)\.(species|genus)_appraise\.tsv$")

frames = []
for path in sorted(glob.glob(os.path.join(results_dir, "*_appraise.tsv"))):
    m = name_re.search(os.path.basename(path))
    if not m:
        print(f"Skipped (name not recognised): {path}")
        continue
    hab = HABITAT[m.group(1)]
    if HAB_NAME[m.group(3)] != hab:
        print(f"WARNING: {os.path.basename(path)} is named '{m.group(3)}' but numbered as {hab}; using {hab}")
    df = pd.read_csv(path, sep="\t", comment="#")
    df = df.assign(Habitat=hab, Reads=READS[m.group(4)], Reference=REF[m.group(5)],
                   Rank=m.group(6).capitalize(), Domain=df["domain"].str.replace("d__", ""))
    frames.append(df)

data   = pd.concat(frames, ignore_index=True)
pooled = data[data["sample"] == "total"]
data   = data[~data["sample"].isin(["total", "average"])]
data.to_csv(out_prefix + ".tsv", sep="\t", index=False)

groups = [(d, r) for d in ["Archaea", "Bacteria"] for r in ["Species", "Genus"]]
offset = {"Sediment": -0.18, "Seawater": 0.18}
rng    = np.random.default_rng(1)

fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharey=True)
for i, reads in enumerate(READS.values()):
    for j, ref in enumerate(REF.values()):
        ax = axes[i, j]
        for k, (dom, rank) in enumerate(groups):
            for hab in ["Sediment", "Seawater"]:
                sel = (data.Reads == reads) & (data.Reference == ref) & (data.Domain == dom) & \
                      (data.Rank == rank) & (data.Habitat == hab)
                y = data.loc[sel, "percent_binned"].to_numpy()
                x = k + offset[hab] + rng.uniform(-0.1, 0.1, len(y))
                ax.scatter(x, y, s=8, alpha=0.5, color=COLOURS[hab], edgecolors="none",
                           label=hab if (k == 0 and i == 0 and j == 0) else None)
                p = pooled[(pooled.Reads == reads) & (pooled.Reference == ref) & (pooled.Domain == dom) &
                           (pooled.Rank == rank) & (pooled.Habitat == hab)]["percent_binned"]
                if len(p):
                    ax.hlines(p.iloc[0], k + offset[hab] - 0.13, k + offset[hab] + 0.13, color="black", lw=2)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([f"{d}\n{r}" for d, r in groups])
        ax.set_title(f"{reads}, against {ref}", fontsize=11)
        ax.set_ylim(-3, 103)
        ax.grid(axis="y", color="0.9")
        ax.set_axisbelow(True)
        if j == 0:
            ax.set_ylabel("mOTU coverage assigned (%)")

handles, labels = axes[0, 0].get_legend_handles_labels()
handles.append(plt.Line2D([], [], color="black", lw=2))
labels.append("Pooled over samples")
fig.legend(handles, labels, loc="lower center", ncol=3, frameon=False, markerscale=2)
fig.tight_layout(rect=(0, 0.05, 1, 1))
for ext in ("png", "svg"):
    fig.savefig(f"{out_prefix}.{ext}", dpi=300)

print(data.groupby(["Reads", "Reference", "Habitat", "Domain", "Rank"])["percent_binned"]
          .agg(["count", "mean", "median"]).round(1).to_string())
