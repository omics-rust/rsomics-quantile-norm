# rsomics-quantile-norm

Quantile-normalize a gene x sample count matrix so every sample ends up with
an identical value distribution — the standard preprocessing step before
microarray/RNA-seq differential analysis (limma `normalizeQuantiles`).

```
rsomics-quantile-norm counts.tsv -o normalized.tsv
```

Input is a TSV: a header row of sample IDs (first cell is the gene-id column
label) and one row per gene (`gene_id<TAB>v1<TAB>v2…`). Output is the
normalized matrix in the same shape.

## Method

For each sample, sort its values; the mean across samples of the values at
each sort position forms the common reference quantile curve. Each value is
then replaced by that curve evaluated at the value's rank. Ties use average
ranks (R `rank(ties.method="average")`) and the reference curve is linearly
interpolated at the resulting fractional rank — so a run of equal values in a
sample maps to the mean of the reference quantiles they collectively span.

## Origin

This crate is an independent Rust reimplementation of limma's
`normalizeQuantiles` based on:

- The published quantile-normalization method (Bolstad, Irizarry, Åstrand &
  Speed, *Bioinformatics* 19(2):185–193, 2003, doi:10.1093/bioinformatics/19.2.185).
- The documented behavior of `limma::normalizeQuantiles` (average-tie ranks,
  uniform-grid linear interpolation of the rank-mean reference).
- Black-box value-level testing against the limma binary (`tests/compat.rs`
  diffs ours against `limma::normalizeQuantiles` and a committed golden).

No source code from the GPL upstream was used as reference during
implementation.

License: MIT OR Apache-2.0.
Upstream credit: [limma](https://bioconductor.org/packages/limma/) (GPL >= 2).
