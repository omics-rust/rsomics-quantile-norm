#!/usr/bin/env Rscript
# Reference oracle: limma::normalizeQuantiles on a gene x sample count matrix.
# Usage: Rscript quantile_norm_oracle.R <counts.tsv> [out.tsv]
# Emits the normalized matrix as TSV (gene id + 6-decimal values), matching
# rsomics-quantile-norm's output formatting for value-level diffing.
suppressMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
counts_path <- args[1]
out_path <- if (length(args) >= 2) args[2] else ""

df <- read.delim(counts_path, header = TRUE, check.names = FALSE,
                 row.names = 1)
genes <- rownames(df)

A <- as.matrix(df)
storage.mode(A) <- "double"
N <- normalizeQuantiles(A, ties = TRUE)

con <- if (nchar(out_path) > 0) file(out_path, "w") else stdout()
# read.delim folds the first header cell into row.names, so echo the raw header.
raw_header <- readLines(counts_path, n = 1)
writeLines(raw_header, con)

vals <- formatC(N, format = "f", digits = 6)
for (i in seq_len(nrow(N))) {
  writeLines(paste(c(genes[i], vals[i, ]), collapse = "\t"), con)
}
if (nchar(out_path) > 0) close(con)
