# CanvasXpressBio

Interactive [CanvasXpress](https://canvasxpress.org) visualizations for core
[Bioconductor](https://bioconductor.org) data structures. `CanvasXpressBio` complements the
CRAN `canvasXpress` package with methods that take a Bioconductor object directly and render
it as a reproducible CanvasXpress widget.

## Install

Once accepted into Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("CanvasXpressBio")
```

## Use

```r
library(CanvasXpressBio)
library(SummarizedExperiment)

# Heatmap from a SummarizedExperiment (colData -> sample annotations,
# rowData -> feature annotations). Works for any class that extends
# SummarizedExperiment, including SingleCellExperiment and DESeqDataSet.
cxplot(se, n = 25, title = "Top 25 variable features")

# Volcano plot from a differential-expression result table.
cxvolcano(DESeq2::results(dds))
```

See the vignette (`vignette("CanvasXpressBio")`) for a worked example.

## Relationship to the CRAN package

The general-purpose R interface to CanvasXpress lives on CRAN as
[`canvasXpress`](https://CRAN.R-project.org/package=canvasXpress). `CanvasXpressBio` depends on
it and adds only the Bioconductor-object methods, so the two are installed side by side.

## License

GPL-3. The underlying CanvasXpress library is dual-licensed; see
<https://canvasxpress.org/license.html>.
