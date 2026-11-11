#!/usr/bin/python3
"""
clustomatic.modified.py

Modified from clustomatic.py of clust-o-matic
(https://github.com/Helmholtz-HIPS/clustomatic_source). The scoring, linkage and
clustering are unchanged. Changes:
  - two extra arguments: an output folder and a number of DIAMOND threads
  - DIAMOND makedb and blastp run with --threads
  - progress messages
  - results written to <outdir>/cluster_assignments.txt (header BGC, GCF) and
    <outdir>/distance_matrix.txt (condensed distance matrix, one value per line)
    instead of printing the clusters to standard output

NOTE: the original modified script was not kept. This version is reconstructed
from the published script and the messages and outputs of the run
("Distance matrix written to ...", "Cluster assignments written to ...").
The distance matrix is written before linkage here, because
fastcluster.linkage(preserve_input=False) overwrites the matrix in place.

Usage:
    python clustomatic.modified.py input.fasta threshold outdir threads
Protein headers must be >clusterID_proteinID with a single underscore.
"""
import scipy
import scipy.cluster
import os
import sys
import numpy as np
import fastcluster
import shutil
from Bio import SeqIO


def parse_input(path, lib):
    index = 0
    for fasta in SeqIO.parse(open(path), "fasta"):
        name, sequence = fasta.id, str(fasta.seq)
        (qcid, qgid) = name.split("_")
        if not lib.get(qcid):
            lib[qcid] = {"index": index}
            index = index + 1
        lib[qcid]["length"] = lib[qcid].get("length", 0) + len(sequence)


def prepare_matrix(lib):
    return np.ones(int(((len(lib) - 1) * (len(lib))) / 2))


def make_blast(path, threads):
    pid = os.posix_spawnp("diamond", [
        "diamond", "makedb", "--quiet",
        "--threads", str(threads),
        "--db", "{0}.db".format(path),
        "--in", "{0}".format(path)
    ], os.environ)
    os.waitpid(pid, 0)


def run_blast(path, threads):
    pid = os.posix_spawnp("diamond", [
        "diamond", "blastp", "--quiet",
        "--threads", str(threads),
        "--db", "{0}.db".format(path),
        "-q", "{0}".format(path),
        "--max-target-seqs", "0",
        "--outfmt", "6", "qseqid", "sseqid", "nident",
        "--out", "{0}.output".format(path)
    ], os.environ)
    os.waitpid(pid, 0)


def parse_blast(path, matrix, lib):
    n = len(lib)
    with open("{0}.output".format(path), "r") as output:
        for pool in generate_query_pools(output):
            track = {}
            q = lib[pool[0]['qcid']]
            i = q['index']
            for line in sorted(pool, key=lambda x: -x['score']):
                if not track.setdefault(line['qcid'] + line['scid'], {}).get(line['qgid']) and \
                   not track.setdefault(line['qcid'] + line['scid'], {}).get(line['sgid']):
                    track[line['qcid'] + line['scid']][line['qgid']] = True
                    track[line['qcid'] + line['scid']][line['sgid']] = True
                    h = lib[line['scid']]
                    j = h['index']
                    if i < j:
                        k = int(n * (n - 1) / 2 - (n - i) * ((n - i) - 1) / 2 + j - i - 1)
                        matrix[k] = matrix[k] - 2 * line['score'] / (q['length'] + h['length'])
                        if matrix[k] < 0:
                            matrix[k] = 0


def generate_query_pools(stream):
    pool = []
    for line in stream:
        line = line.split()
        (qcid, qgid) = line[0].split('_')
        (scid, sgid) = line[1].split('_')
        if len(pool) > 0 and pool[-1]['qcid'] != qcid:
            yield pool
            pool = []
        pool.append({'qcid': qcid, 'qgid': qgid, 'scid': scid, 'sgid': sgid, 'score': int(line[2])})
    yield pool


def build_linkage(matrix):
    linkage = fastcluster.linkage(matrix, method='average', preserve_input=False)
    assert scipy.cluster.hierarchy.is_valid_linkage(linkage), 'reconstructed linkage not valid'
    return linkage


def form_clusters(linkage, threshold):
    return scipy.cluster.hierarchy.fcluster(linkage, t=threshold, criterion="distance")


def clustomatic(path, threshold, outdir, threads):
    os.makedirs(outdir, exist_ok=True)
    lib = {}

    print("RUNNING CLUSTOMATIC", flush=True)
    parse_input(path, lib)
    make_blast(path, threads)

    print("RUNNING BLAST", flush=True)
    run_blast(path, threads)

    print("PARSING BLAST", flush=True)
    matrix = prepare_matrix(lib)
    parse_blast(path, matrix, lib)

    dist_out = os.path.join(outdir, "distance_matrix.txt")
    np.savetxt(dist_out, matrix, fmt="%.6f")   # before linkage, which overwrites matrix

    print("GETTING LINKAGE", flush=True)
    linkage = build_linkage(matrix)

    print("GETTING CLUSTERS", flush=True)
    clusters = form_clusters(linkage, threshold)

    print("WRITING OUTPUT", flush=True)
    clust_out = os.path.join(outdir, "cluster_assignments.txt")
    with open(clust_out, "w") as out:
        out.write("BGC\tGCF\n")
        for (i, key) in enumerate(lib):
            out.write("{0}\t{1}\n".format(key, clusters[i]))

    print("Distance matrix written to {0}".format(dist_out))
    print("Cluster assignments written to {0}".format(clust_out))


def main():
    if len(sys.argv) < 5:
        raise SystemExit("Usage: clustomatic.modified.py input.fasta threshold outdir threads\n"
                         "Fasta protein records must have following naming: >clusterID_proteinID")
    path, threshold, outdir, threads = sys.argv[1], float(sys.argv[2]), sys.argv[3], int(sys.argv[4])
    if not os.path.isfile(path):
        raise SystemExit("File not found: {0}".format(path))
    if threshold < 0 or threshold > 1:
        raise SystemExit("Threshold has to be a float between 0 and 1")
    if not shutil.which("diamond"):
        raise SystemExit("Please install diamond (https://github.com/bbuchfink/diamond)")
    clustomatic(path, threshold, outdir, threads)


if __name__ == "__main__":
    main()
