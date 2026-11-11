#!/usr/bin/env python3
"""
5.07_enrichment_mannwhitney.py

MAG annotation, step 7: compare gene copy numbers between seawater and
sediment MAGs (796 highly complete MAGs: 449 seawater, 347 sediment), one
annotation at a time, with two-sided Mann-Whitney U tests and
Benjamini-Hochberg correction.


Input
    --annotation  long table from 5.06 (Genome <tab> Annotation; one row per gene)
    --groups      MAGs to compare (columns Entry and Source; Source is
                  Seawater or Sediment)
Output
    --out         one row per annotation: group means, log2 fold changes,
                  U statistic, p-value, and BH-adjusted p-value

Usage
    python 2_genome_centric/5_MAGsAnnot/scripts/5.07_enrichment_mannwhitney.py \
        --annotation annotations/ko.tsv --groups genome_group.tsv --out ko_enrichment.csv
"""

import argparse

import numpy as np
import pandas as pd
import scipy.stats as stats

G1, G2 = "Seawater", "Sediment"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--groups", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    ann = pd.read_csv(a.annotation, sep="\t")
    groups = pd.read_csv(a.groups, sep="\t")
    g1 = groups.loc[groups["Source"] == G1, "Entry"].tolist()
    g2 = groups.loc[groups["Source"] == G2, "Entry"].tolist()
    print(f"{G1}: {len(g1)} MAGs; {G2}: {len(g2)} MAGs")

    # Count matrix (annotations x MAGs)
    counts = (ann[ann["Genome"].isin(g1 + g2)]
              .groupby(["Annotation", "Genome"]).size()
              .unstack(fill_value=0))
    missing = set(g1 + g2) - set(counts.columns)
    if missing:
        print(f"Warning: {len(missing)} MAGs have no annotation of this type and are not tested")
    g1 = [m for m in g1 if m in counts.columns]
    g2 = [m for m in g2 if m in counts.columns]
    print(f"{counts.shape[0]} annotations")

    rows = []
    for gene, row in counts.iterrows():
        x1, x2 = row[g1].to_numpy(), row[g2].to_numpy()
        m1, m2 = x1.mean(), x2.mean()
        u, p = stats.mannwhitneyu(x1, x2, alternative="two-sided")
        fc = np.log2((m1 + 1e-99) / (m2 + 1e-99))
        rows.append({"annotation": gene,
                     "group_1": G1, "group_2": G2,
                     "enriched_in": G1 if m1 > m2 else (G2 if m2 > m1 else "NA"),
                     "group_1_mean": m1, "group_2_mean": m2,
                     "log2FC_g1/g2": fc, "log2FC_g2/g1": -fc,
                     "U": u, "pvalue": p})

    res = pd.DataFrame(rows)
    res["corrected_pvalue"] = stats.false_discovery_control(res["pvalue"], method="bh")
    res.sort_values("corrected_pvalue").to_csv(a.out, index=False)
    print(f"Adjusted p <= 0.05: {(res['corrected_pvalue'] <= 0.05).sum()}")


if __name__ == "__main__":
    main()
