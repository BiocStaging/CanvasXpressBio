# CanvasXpressBio 0.99.0

* Initial submission to Bioconductor.
* `cxplot()` generic renders a `SummarizedExperiment` (and extending classes such as
  `SingleCellExperiment` and `DESeqDataSet`) as an interactive CanvasXpress heatmap, with
  `colData` as sample annotations and `rowData` as feature annotations. `n=` keeps only the
  most variable features.
* `cxvolcano()` draws a volcano plot from a differential-expression result table
  (`DESeqResults`, `limma::topTable()`, or a data frame).
