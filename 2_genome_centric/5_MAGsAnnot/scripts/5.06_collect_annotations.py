#!/usr/bin/env python3
"""
5.06_collect_annotations.py

MAG annotation, step 6: collect per-MAG annotations (5.01 to 5.04) into one
long table per annotation type, with one row per annotated gene:
    Genome <tab> Annotation

Annotation types written to --outdir:
    ko.tsv            KOfam (DRAM ko_id)
    pfam.tsv          Pfam (DRAM pfam_hits; one row per Pfam domain)
    merops.tsv        MEROPS peptidases (DRAM peptidase_id)
    cazy.tsv          CAZyme families (run_dbcan overview.tsv, recommended
                      family; one row per family for multi-family proteins)
    substrate.tsv     CGC substrates from dbCAN-PUL (run_dbcan
                      substrate_prediction.tsv, column 3; one row per CGC)
    sulfatlas.tsv     SulfAtlas families (from the DIAMOND subject ID)
    hyddb.tsv         HydDB hydrogenase subgroups (filtered hits)
    marker_genes.tsv  44 metabolic marker genes (best hit per protein)

Usage:
    python 2_genome_centric/5_MAGsAnnot/scripts/5.06_collect_annotations.py \
        --dram   <5_MAG_annotation/1_DRAM> \
        --dbcan  <5_MAG_annotation/2_CAZymes_CGCs> \
        --sulf   <5_MAG_annotation/3_sulphatases> \
        --hyddb  <5_MAG_annotation/4_hydrogenases_marker_genes/hyddb> \
        --marker <5_MAG_annotation/4_hydrogenases_marker_genes/marker_best_hits> \
        --outdir <5_MAG_annotation/6_enrichment_tests/annotations>
"""

import argparse
from pathlib import Path

import pandas as pd


def write(df, path):
    df = df.dropna(subset=["Annotation"])
    df = df[df["Annotation"].astype(str).str.strip() != ""]
    df.to_csv(path, sep="\t", index=False)
    print(f"{path.name}: {len(df)} genes in {df['Genome'].nunique()} genomes")


def collect_dram(dram_dir, outdir):
    """KO, Pfam, and MEROPS from DRAM annotations.tsv (one folder per MAG)."""
    ko, pfam, merops = [], [], []
    for d in sorted(Path(dram_dir).glob("*__annotation")):
        genome = d.name.removesuffix("__annotation")
        ann = pd.read_csv(d / "annotations.tsv", sep="\t", index_col=0)
        if "ko_id" in ann:
            ko.append(pd.DataFrame({"Genome": genome, "Annotation": ann["ko_id"]}))
        if "pfam_hits" in ann:
            p = ann["pfam_hits"].dropna().str.split("; ").explode()
            pfam.append(pd.DataFrame({"Genome": genome, "Annotation": p}))
        if "peptidase_id" in ann:
            merops.append(pd.DataFrame({"Genome": genome, "Annotation": ann["peptidase_id"]}))
    write(pd.concat(ko), outdir / "ko.tsv")
    write(pd.concat(pfam), outdir / "pfam.tsv")
    write(pd.concat(merops), outdir / "merops.tsv")


def collect_dbcan(dbcan_dir, outdir):
    """CAZyme families: recommended result (column 7) of run_dbcan overview.tsv."""
    rows = []
    for d in sorted(p for p in Path(dbcan_dir).iterdir() if p.is_dir()):
        genome = d.name.removesuffix("__annotation")
        ov = pd.read_csv(d / "overview.tsv", sep="\t")
        fam = ov.iloc[:, 6]
        fam = fam[fam != "-"].str.split("|").explode().str.replace(r"\.hmm$", "", regex=True)
        rows.append(pd.DataFrame({"Genome": genome, "Annotation": fam}))
    write(pd.concat(rows), outdir / "cazy.tsv")

    # dbCAN-PUL substrate of each CGC (column 3 of substrate_prediction.tsv)
    rows = []
    for d in sorted(p for p in Path(dbcan_dir).iterdir() if p.is_dir()):
        f = d / "substrate_prediction.tsv"
        if not f.exists():
            continue
        sub = pd.read_csv(f, sep="\t")
        if sub.shape[0] == 0:
            continue
        rows.append(pd.DataFrame({"Genome": d.name.removesuffix("__annotation"), "Annotation": sub.iloc[:, 2]}))
    write(pd.concat(rows), outdir / "substrate.tsv")


def collect_sulfatlas(sulf_dir, outdir):
    """SulfAtlas family from the DIAMOND subject ID (TODO: confirm ID format)."""
    rows = []
    for f in sorted(Path(sulf_dir).glob("*.tsv")):
        genome = f.stem.removesuffix("__annotation")
        if f.stat().st_size == 0:
            continue
        hits = pd.read_csv(f, sep="\t", header=None)
        family = hits[1].astype(str).str.split("|").str[-1]
        rows.append(pd.DataFrame({"Genome": genome, "Annotation": family}))
    write(pd.concat(rows), outdir / "sulfatlas.tsv")


def collect_hyddb(hyddb_dir, outdir):
    rows = []
    for f in sorted(Path(hyddb_dir).glob("*.tsv")):
        hits = pd.read_csv(f, sep="\t")
        rows.append(pd.DataFrame({"Genome": f.stem, "Annotation": hits["closest_hydrogenase_group"]}))
    write(pd.concat(rows), outdir / "hyddb.tsv")


def collect_markers(marker_dir, outdir):
    """Marker gene database name: last column of the best-hit tables."""
    rows = []
    for f in sorted(Path(marker_dir).glob("*.tsv")):
        if f.stat().st_size == 0:
            continue
        hits = pd.read_csv(f, sep="\t", header=None)
        rows.append(pd.DataFrame({"Genome": f.stem, "Annotation": hits[hits.columns[-1]]}))
    write(pd.concat(rows), outdir / "marker_genes.tsv")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dram", required=True)
    ap.add_argument("--dbcan", required=True)
    ap.add_argument("--sulf", required=True)
    ap.add_argument("--hyddb", required=True)
    ap.add_argument("--marker", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()

    outdir = Path(a.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    collect_dram(a.dram, outdir)
    collect_dbcan(a.dbcan, outdir)
    collect_sulfatlas(a.sulf, outdir)
    collect_hyddb(a.hyddb, outdir)
    collect_markers(a.marker, outdir)


if __name__ == "__main__":
    main()
