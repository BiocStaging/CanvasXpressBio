# CanvasXpressBio 0.99.3

* Version bump to re-run the Bioconductor staging build (no code changes; 0.99.2 built
  and passed checks but failed to publish).

# CanvasXpressBio 0.99.2

* Vignette charts now render in the built HTML: require `canvasXpress (>= 1.70.3)`, which
  ships a CanvasXpress library that BiocStyle's self-contained HTML post-processing no longer
  corrupts (earlier versions were mangled on inlining, leaving the widgets blank). Verified by
  rendering the vignette against the CRAN 1.70.3 library.

# CanvasXpressBio 0.99.1

* Add maintainer ORCID iD to `Authors@R` (BiocCheck note).

# CanvasXpressBio 0.99.0

* Initial submission to Bioconductor.
* `cxplot()` generic renders a `SummarizedExperiment` (and extending classes such as
  `SingleCellExperiment` and `DESeqDataSet`) as an interactive CanvasXpress heatmap, with
  `colData` as sample annotations and `rowData` as feature annotations. `n=` keeps only the
  most variable features.
* `cxvolcano()` draws a volcano plot from a differential-expression result table
  (`DESeqResults`, `limma::topTable()`, or a data frame).
