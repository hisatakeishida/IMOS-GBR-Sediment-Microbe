#!/usr/bin/env python3
"""
Combine per-sample SingleM taxonomic profiles into one coverage matrix.

Input : a directory of SingleM `pipe -p` profiles (<sample>_prof.tsv), each with
        the columns: sample, coverage, taxonomy
Output: CSV with one row per sample and one column per taxon (coverage);
        taxa absent from a sample are 0

Usage:
    python 1.3_singlem_prof_coverage.py \
        --indir <path/to/prof> \
        --out singlem_10Mseqs.profsum.RS226.csv
"""

import argparse
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--indir", required=True, help="directory of *_prof.tsv files")
    parser.add_argument("--out", required=True, help="output CSV path")
    parser.add_argument("--pattern", default="*_prof.tsv", help="profile file pattern")
    args = parser.parse_args()

    files = sorted(Path(args.indir).glob(args.pattern))
    if not files:
        raise SystemExit(f"No files matching {args.pattern} in {args.indir}")

    profiles, samples = [], []
    for f in files:
        df = pd.read_csv(f, sep="\t", usecols=["sample", "coverage", "taxonomy"])
        if df.empty:
            # Keep samples with no detected taxa as an all-zero row
            name = f.name.removesuffix("_prof.tsv")
            print(f"Warning: {f.name} has no taxa; kept as all-zero row ({name})")
        else:
            name = df["sample"].iloc[0]
            profiles.append(df)
        samples.append(name)

    long = pd.concat(profiles, ignore_index=True)
    matrix = (long.pivot_table(index="sample", columns="taxonomy", values="coverage",
                               aggfunc="sum", fill_value=0)
                  .reindex(sorted(samples), fill_value=0)
                  .sort_index(axis=1))
    matrix.index.name = "sample"

    matrix.to_csv(args.out)
    print(f"{matrix.shape[0]} samples x {matrix.shape[1]} taxa written to {args.out}")


if __name__ == "__main__":
    main()