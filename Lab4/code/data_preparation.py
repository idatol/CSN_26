"""
For every language in list.txt, computes the metrics of the syntactic
dependency tree of every sentence in the UD treebanks and writes them to
data/<language>_dependency_tree_metrics.txt.

File format (no header, space separated), one row per sentence:
    n  k2  d
where
    n  = number of vertices (tokens) of the tree
    k2 = <k^2>, degree 2nd moment about zero = (1/n) * sum_i k_i^2
    d  = <d>, mean edge length = (1/(n-1)) * sum_edges |pos(dependent) - pos(head)|

Decisions:
  - All treebanks of a language are used (same as Lab 2 and 3), e.g. all
    UD_English-* folders and all their .conllu files (train/dev/test).
  - Multiword tokens (ids like 3-4) and empty nodes (ids like 5.1) are skipped,
    so a sentence keeps its regular tokens 1..n in their linear order.
  - Punctuation is kept (as in previous labs), so the tree stays connected.
  - The tree is treated as undirected: k_i counts the head and the dependents
    of word i. The edge from the root to the artificial node 0 is not an edge.
  - Sentences with n < 2 are skipped because <d> is undefined (no edges).
  - Sentences that do not form a valid tree (ids not 1..n, a head outside the
    sentence, not exactly one root, or a cycle) are skipped and counted.

Run from the Lab4/code folder:  python3 data_preparation.py
"""
import glob
import os

UD_DIR = "../ud-treebanks-v2.18"
DATA_DIR = "../data"
os.makedirs(DATA_DIR, exist_ok=True)

with open("../list.txt", encoding="utf-8") as f:
    languages = [line.split()[0] for line in f if line.strip() and not line.startswith("language")]
print(f"Languages: {', '.join(languages)}")


def read_sentences(path):
    """Yields the list of heads of every sentence: heads[i-1] is the head of token i."""
    ids, heads = [], []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:  # blank line ends a sentence
                if ids:
                    yield ids, heads
                ids, heads = [], []
                continue
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if "-" in cols[0] or "." in cols[0]:  # multiword token or empty node
                continue
            ids.append(int(cols[0]))
            heads.append(int(cols[6]) if cols[6].isdigit() else -1)
    if ids:
        yield ids, heads


def is_tree(ids, heads):
    n = len(ids)
    if ids != list(range(1, n + 1)):
        return False
    if any(h < 0 or h > n for h in heads) or heads.count(0) != 1:
        return False
    # every vertex must reach the root without repeating a vertex (no cycles)
    for start in range(1, n + 1):
        v, steps = start, 0
        while v != 0:
            v = heads[v - 1]
            steps += 1
            if steps > n:
                return False
    return True


def tree_metrics(heads):
    n = len(heads)
    degree = [0] * (n + 1)
    total_length = 0
    for dep, head in enumerate(heads, start=1):
        if head == 0:
            continue
        degree[dep] += 1
        degree[head] += 1
        total_length += abs(dep - head)
    k2 = sum(k * k for k in degree[1:]) / n
    d = total_length / (n - 1)
    q = sum(k % 2 for k in degree[1:])  # vertices of odd degree, needed for bound (2)
    return n, k2, d, q


def satisfies_bounds(n, k2, d, q, eps=1e-9):
    """Bounds (1) and (2) of the lab guide, checked with the exact q."""
    ok1 = 4 - 6 / n - eps <= k2 <= n - 1 + eps
    ok2 = (n * k2 + q) / 4 - eps <= (n - 1) * d <= (3 * (n - 1) ** 2 + 1 - n % 2) / 4 + eps
    return ok1 and ok2


for language in languages:
    files = sorted(glob.glob(os.path.join(UD_DIR, f"UD_{language}-*", "*.conllu")))
    rows, too_short, invalid, out_of_bounds = [], 0, 0, 0
    for path in files:
        for ids, heads in read_sentences(path):
            if len(ids) < 2:
                too_short += 1
            elif not is_tree(ids, heads):
                invalid += 1
            else:
                n, k2, d, q = tree_metrics(heads)
                out_of_bounds += not satisfies_bounds(n, k2, d, q)
                rows.append((n, k2, d))

    out = os.path.join(DATA_DIR, f"{language}_dependency_tree_metrics.txt")
    with open(out, "w", encoding="utf-8") as f:
        for n, k2, d in rows:
            f.write(f"{n} {k2!r} {d!r}\n")
    print(f"{language}: {len(files)} files, {len(rows)} trees "
          f"(skipped {too_short} with n < 2, {invalid} invalid; "
          f"{out_of_bounds} violate the bounds) -> {out}")
