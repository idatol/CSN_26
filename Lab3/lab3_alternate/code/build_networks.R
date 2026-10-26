# Builds the undirected global syntactic dependency networks from the UD treebanks.
# Usage: Rscript build_networks.R <path/to/ud-treebanks-v2.18> <output_dir>
#
# Vertices are lowercased word forms. There is an edge {head, dependent} if the pair
# occurs in any treebank of the language, regardless of direction. Punctuation
# (UPOS = PUNCT) is ignored. Loops and repeated edges are removed.
# Output per language: <language>_dependency_network.txt
#   line 1: "N E"; then one edge per line: word1<TAB>word2
#   (a tab, because some word forms contain spaces).
# Also writes network_summary.txt (Table 1: language N E <k> delta).

library(data.table)

# byte-wise string comparison: locale collation can treat different strings as equal,
# which would turn an edge into a loop in the canonical ordering below
Sys.setlocale("LC_COLLATE", "C")

args <- commandArgs(trailingOnly = TRUE)
ud_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

languages <- c("Arabic", "Basque", "Catalan", "Chinese", "Czech",
               "English", "Greek", "Hungarian", "Italian", "Turkish")

# Excluded: no word forms (NYUAD, GUMReddit), historical (Old Italian),
# same text as Chinese-GSD (GSDSimp)
excluded <- c("UD_Arabic-NYUAD", "UD_English-GUMReddit", "UD_Italian-Old", "UD_Chinese-GSDSimp")

# (head, dependent) pairs of one .conllu file
read_edges <- function(file) {
  lines <- readLines(file, encoding = "UTF-8")
  lines <- lines[!grepl("^#", lines) & nzchar(lines)]
  cols <- tstrsplit(lines, "\t", fixed = TRUE, keep = c(1, 2, 4, 7))
  d <- data.table(id = cols[[1]], form = tolower(cols[[2]]),
                  upos = cols[[3]], head = cols[[4]])
  d <- d[grepl("^[0-9]+$", id)]                # drop multiword tokens and empty nodes
  d[, sent := cumsum(id == "1")]               # new sentence whenever id == 1
  d[, id := as.integer(id)]
  d[, head := as.integer(head)]
  d[, head_form := form[match(head, id)], by = sent]
  d[, head_upos := upos[match(head, id)], by = sent]
  d[head > 0 & !is.na(head_form) & upos != "PUNCT" & head_upos != "PUNCT",
    .(head = head_form, dep = form)]
}

summary <- list()

for (lang in languages) {
  treebanks <- setdiff(basename(Sys.glob(file.path(ud_dir, paste0("UD_", lang, "-*")))), excluded)
  files <- unlist(lapply(treebanks, function(t) Sys.glob(file.path(ud_dir, t, "*.conllu"))))
  edges <- rbindlist(lapply(files, read_edges))
  edges <- edges[head != dep]                  # remove loops
  # undirected: store each pair once, in canonical order
  edges <- unique(edges[, .(a = pmin(head, dep), b = pmax(head, dep))])
  N <- length(unique(c(edges$a, edges$b)))
  E <- nrow(edges)
  writeLines(c(paste(N, E), paste(edges$a, edges$b, sep = "\t")),
             file.path(out_dir, paste0(lang, "_dependency_network.txt")))
  summary[[lang]] <- data.frame(language = lang, N = N, E = E,
                                k = 2 * E / N, delta = 2 * E / (N * (N - 1)))
  cat(lang, ":", length(treebanks), "treebanks, N =", N, " E =", E, "\n")
}

write.table(do.call(rbind, summary), file.path(out_dir, "network_summary.txt"),
            quote = FALSE, row.names = FALSE, sep = "\t")
