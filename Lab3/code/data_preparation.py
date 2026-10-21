"""
Builds the undirected global syntactic dependency network of every language in
list.txt and writes it to data/<language>_dependency_network.txt.

File format:
    line 1:   N E
    next E lines: word1<TAB>word2      (one undirected edge per line)

Tab is used as separator because it can never appear inside a CoNLL-U field,
whereas a few word forms contain spaces.

Decisions (same as Lab 2): word forms, lowercased, punctuation kept,
multiword tokens and empty nodes skipped. Loops are removed. A word that only
took part in loops stays as an isolated vertex (counted in N, degree 0).

Run from the Lab3/code folder:  python3 data_preparation.py
"""
import glob
import os

from conllu import parse_incr

UD_DIR = "../ud-treebanks-v2.18"
DATA_DIR = "../data"
os.makedirs(DATA_DIR, exist_ok=True)

with open("../list.txt", encoding="utf-8") as f:
    languages = [line.split()[0] for line in f if line.strip() and not line.startswith("language")]
print(f"Languages: {', '.join(languages)}")


def build_network(paths, use="form", lower=True):
    vertices = set()
    edges = set()  # undirected edges stored as sorted tuples, so no multiedges
    for path in paths:
        with open(path, encoding="utf-8") as f:
            for sentence in parse_incr(f):
                # keep only real tokens: skip multiword tokens and empty nodes
                tokens = {t["id"]: t for t in sentence if isinstance(t["id"], int)}
                for tok in tokens.values():
                    head = tokens.get(tok["head"])
                    if head is None:  # root (head 0) or missing head
                        continue
                    h, d = head[use], tok[use]
                    if h is None or d is None or h == "_" or d == "_":
                        continue
                    if lower:
                        h, d = h.lower(), d.lower()
                    vertices.add(h)
                    vertices.add(d)
                    if h != d:  # remove loops
                        edges.add((h, d) if h < d else (d, h))
    return vertices, edges


rows = []
for language in languages:
    files = sorted(glob.glob(os.path.join(UD_DIR, f"UD_{language}-*", "*.conllu")))
    print(f"{language}: {len(files)} files")
    vertices, edges = build_network(files)

    out = os.path.join(DATA_DIR, f"{language}_dependency_network.txt")
    with open(out, "w", encoding="utf-8") as f:
        f.write(f"{len(vertices)} {len(edges)}\n")
        for u, v in sorted(edges):
            f.write(f"{u}\t{v}\n")

    N, E = len(vertices), len(edges)
    rows.append((language, N, E, 2 * E / N, 2 * E / (N * (N - 1))))
    print(f"  N={N} E={E} -> {out}")

# Table 1
os.makedirs("../results", exist_ok=True)
with open("../results/table1.tsv", "w") as f:
    f.write("language\tN\tE\tmean_k\tdensity\n")
    for lang, N, E, k, delta in rows:
        f.write(f"{lang}\t{N}\t{E}\t{k:.4f}\t{delta:.3e}\n")
print("Table 1 written to ../results/table1.tsv")
