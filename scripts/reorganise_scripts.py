#!/usr/bin/env python3
"""
reorganise_scripts.py

Move the numbered scripts from the top-level scripts/ folder into a scripts/
subfolder of each section directory, and update every reference to them:
  - markdown links in the pages (relative to each page)
  - paths in run commands and script headers (relative to the repository root,
    e.g. "sbatch scripts/1.01_x.sh" becomes "sbatch 1_gene_centric/1_SingleM/scripts/1.01_x.sh")
  - "<path/to/scripts/x>" placeholders inside scripts
It also repairs links between pages that no longer resolve (for example
../4_MAG_recovery/4_MAG_recovery.md), by matching the page file name, or the
page number (5_MAG_annotation.md -> 5_MAG_annotate.md). Files under archive/
are not touched.

Run from the repository root:
    python3 reorganise_scripts.py            # dry run: lists what would change
    python3 reorganise_scripts.py --apply    # move files (git mv) and edit pages
"""
import os
import re
import shutil
import subprocess
import sys

APPLY = "--apply" in sys.argv
ROOT = os.getcwd()
OLD_DIR = "scripts"

# Section number -> section directory (scripts go into <dir>/scripts/)
SECTIONS = {
    "1": "1_gene_centric/1_SingleM",
    "2": "1_gene_centric/2_Protein",
    "3": "1_gene_centric/3_Sylph",
    "4": "2_genome_centric/4_MAGs",
    "5": "2_genome_centric/5_MAGsAnnot",
    "6": "3_combine/6_Mixomics",
    "7": "3_combine/7_Gene_Genome_analysis",
}
# Unnumbered helper scripts
HELPERS = {
    "clustomatic.modified.py": "5",
    "perform_l2norm_clustering.py": "5",
}
SKIP_DIRS = {"archive", ".git"}


def section_of(name):
    if name in HELPERS:
        return HELPERS[name]
    m = re.match(r"(\d+)\.\d+_", name)
    return m.group(1) if m else None


def walk(exts):
    for d, dirs, files in os.walk(ROOT):
        dirs[:] = [x for x in dirs if x not in SKIP_DIRS]
        for f in files:
            if f.endswith(exts):
                yield os.path.join(d, f)


# ---- 1. Plan the moves -----------------------------------------------------------------
old_dir = os.path.join(ROOT, OLD_DIR)
if not os.path.isdir(old_dir):
    sys.exit(f"No {OLD_DIR}/ folder here; run from the repository root.")

new_path = {}   # script name -> new path relative to ROOT
for name in sorted(os.listdir(old_dir)):
    sec = section_of(name)
    if sec is None or sec not in SECTIONS:
        print(f"WARNING: not moved (no section): {OLD_DIR}/{name}")
        continue
    new_path[name] = f"{SECTIONS[sec]}/scripts/{name}"

if os.path.exists(os.path.join(old_dir, "4.09_gtdbtk_classify.sh")):
    print("NOTE: 4.09_gtdbtk_classify.sh is the old version; replace it with "
          "4.09_gtdbtk_taxonomy_phylogeny.sh (the name used in 4_MAG_recovery.md).")

# ---- 2. Rewrite references -------------------------------------------------------------
pages = {}  # basename -> path relative to ROOT
for p in walk((".md",)):
    pages.setdefault(os.path.basename(p), os.path.relpath(p, ROOT))


def resolve_page(base):
    if base in pages:
        return pages[base]
    m = re.match(r"(\d+)_", base)
    if m:
        hits = [v for k, v in pages.items() if k.startswith(m.group(1) + "_")
                and os.path.dirname(v) in SECTIONS.values()]
        if len(hits) == 1:
            return hits[0]
    return None


names_re = "|".join(re.escape(n) for n in sorted(new_path, key=len, reverse=True))
plain_re = re.compile(r"(?<!\[`)(?<![\w./<-])" + OLD_DIR + r"/(" + names_re + r")\b")
placeholder_re = re.compile(r"<path/to/" + OLD_DIR + r"/(" + names_re + r")>")
link_re = re.compile(r"\]\(([^)\s]+)\)")


def fix_links(text, md_rel):
    md_dir = os.path.dirname(md_rel)
    sec_dir = md_dir if md_dir in SECTIONS.values() else None

    def repl(m):
        target = m.group(1)
        if re.match(r"[a-z]+://|#|mailto:", target):
            return m.group(0)
        path, _, anchor = target.partition("#")
        anchor = "#" + anchor if anchor else ""
        resolved = os.path.normpath(os.path.join(md_dir, path))
        new = None
        if resolved == OLD_DIR:                                  # link to the scripts folder
            if sec_dir:
                new = os.path.relpath(f"{sec_dir}/scripts", md_dir) + "/"
        elif resolved.startswith(OLD_DIR + "/"):
            name = resolved[len(OLD_DIR) + 1:]
            if name in new_path:
                new = os.path.relpath(new_path[name], md_dir)
        elif path.endswith(".md") and not os.path.exists(os.path.join(ROOT, resolved)):
            hit = resolve_page(os.path.basename(path))
            if hit:
                new = os.path.relpath(hit, md_dir)
            else:
                print(f"WARNING: unresolved link in {md_rel}: {target}")
        return m.group(0) if new is None else f"]({new}{anchor})"

    return link_re.sub(repl, text)


def fix_text(text):
    text = placeholder_re.sub(lambda m: f"<path/to/repo>/{new_path[m.group(1)]}", text)
    return plain_re.sub(lambda m: new_path[m.group(1)], text)


edits = {}
for p in walk((".md", ".sh", ".R", ".py")):
    rel = os.path.relpath(p, ROOT)
    if rel == "reorganise_scripts.py":
        continue
    with open(p, encoding="utf-8") as fh:
        old = fh.read()
    new = fix_text(old)
    if rel.endswith(".md"):
        new = new.replace("the top-level [`scripts/`]", "the [`scripts/`]")
        new = fix_links(new, rel)
    if new != old:
        edits[p] = new

# ---- 3. Report or apply ----------------------------------------------------------------
print(f"\n{len(new_path)} scripts to move; {len(edits)} files with updated references\n")
for name, dest in new_path.items():
    print(f"  {OLD_DIR}/{name} -> {dest}")
for p in edits:
    print(f"  edit: {os.path.relpath(p, ROOT)}")

if not APPLY:
    print("\nDry run only. Re-run with --apply to make the changes.")
    sys.exit(0)

for p, text in edits.items():
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(text)

in_git = subprocess.run(["git", "rev-parse"], cwd=ROOT, capture_output=True).returncode == 0
for name, dest in new_path.items():
    os.makedirs(os.path.dirname(os.path.join(ROOT, dest)), exist_ok=True)
    src = f"{OLD_DIR}/{name}"
    if in_git and subprocess.run(["git", "ls-files", "--error-unmatch", src], cwd=ROOT,
                                 capture_output=True).returncode == 0:
        subprocess.run(["git", "mv", src, dest], cwd=ROOT, check=True)
    else:
        shutil.move(os.path.join(ROOT, src), os.path.join(ROOT, dest))

if not os.listdir(old_dir):
    os.rmdir(old_dir)
print("\nDone. Check with: git status && git diff")