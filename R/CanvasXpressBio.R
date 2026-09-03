#' CanvasXpress visualizations for Bioconductor data structures
#'
#' Bridges the CanvasXpress interactive visualization library to core
#' Bioconductor classes. The [cxplot()] generic renders a
#' [SummarizedExperiment::SummarizedExperiment] (and any class that extends it,
#' such as `SingleCellExperiment` or `DESeqDataSet`) as an interactive
#' CanvasXpress heatmap carrying its row and column annotations, and
#' [cxvolcano()] draws a volcano plot from a differential-expression result
#' table.
#'
#' @name CanvasXpressBio-package
#' @keywords internal
#' @importFrom methods setGeneric setMethod
#' @importFrom stats var
#' @importFrom canvasXpress canvasXpress
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
"_PACKAGE"


#' Visualize a Bioconductor object with CanvasXpress
#'
#' Renders a Bioconductor data object as an interactive CanvasXpress chart. The
#' assay matrix becomes the chart data, with `colData` supplied as sample
#' (column) annotations and `rowData` as variable (row) annotations, so the
#' annotations are available for interactive grouping, coloring and sorting.
#'
#' @param object A Bioconductor data object. Currently a
#'   [SummarizedExperiment::SummarizedExperiment] or any class extending it
#'   (e.g. `SingleCellExperiment`, `DESeqDataSet`).
#' @param ... Additional CanvasXpress configuration passed through to
#'   [canvasXpress::canvasXpress()] (for example `title` or `colorScheme`).
#'
#' @return A `canvasXpress` htmlwidget.
#'
#' @examples
#' se <- SummarizedExperiment::SummarizedExperiment(
#'     assays  = list(counts = matrix(rnorm(20), nrow = 5,
#'                    dimnames = list(paste0("g", 1:5), paste0("s", 1:4)))),
#'     colData = data.frame(group = rep(c("A", "B"), 2),
#'                          row.names = paste0("s", 1:4)))
#' cxplot(se, title = "Example heatmap")
#'
#' @export
setGeneric("cxplot", function(object, ...) standardGeneric("cxplot"))


#' @rdname cxplot
#'
#' @param assay Assay to plot, given as a name or an index. Defaults to the
#'   first assay.
#' @param graphType CanvasXpress graph type. Defaults to `"Heatmap"`.
#' @param n Optional integer. When supplied, only the `n` most variable
#'   features (rows) are kept before plotting — a common convenience for large
#'   genomics matrices. Defaults to `NULL`, which keeps all features.
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @export
setMethod(
    "cxplot", "SummarizedExperiment",
    function(object, assay = 1L, graphType = "Heatmap", n = NULL, ...) {
        mat <- as.matrix(SummarizedExperiment::assay(object, assay))
        if (is.null(rownames(mat))) {
            rownames(mat) <- paste0("V", seq_len(nrow(mat)))
        }
        if (is.null(colnames(mat))) {
            colnames(mat) <- paste0("S", seq_len(ncol(mat)))
        }

        var_annot <- as.data.frame(SummarizedExperiment::rowData(object))
        smp_annot <- as.data.frame(SummarizedExperiment::colData(object))
        rownames(smp_annot) <- colnames(mat)
        if (nrow(var_annot) == nrow(mat)) {
            rownames(var_annot) <- rownames(mat)
        }

        # Keep the n most variable features when requested.
        if (!is.null(n) && is.finite(n) && n < nrow(mat)) {
            vars_by_row <- apply(mat, 1L, stats::var)
            keep <- sort(order(vars_by_row, decreasing = TRUE)[seq_len(n)])
            if (nrow(var_annot) == length(vars_by_row)) {
                var_annot <- var_annot[keep, , drop = FALSE]
            }
            mat <- mat[keep, , drop = FALSE]
        }

        use_var_annot <- nrow(var_annot) == nrow(mat) && ncol(var_annot)
        canvasXpress::canvasXpress(
            data      = mat,
            smpAnnot  = if (ncol(smp_annot)) smp_annot else NULL,
            varAnnot  = if (use_var_annot) var_annot else NULL,
            graphType = graphType,
            ...)
    })


#' Volcano plot of differential-expression results with CanvasXpress
#'
#' Draws an interactive volcano plot (log fold change against
#' \eqn{-\log_{10}} p-value) from a table of differential-expression results.
#' The input is coerced with [as.data.frame()], so a `DESeqResults` object, a
#' `limma::topTable()` data frame, or any comparable table works.
#'
#' @param x A data-frame-like table of results.
#' @param logfc,pval Column names holding the log fold change and the (adjusted)
#'   p-value. The defaults suit `DESeq2` output.
#' @param ... Additional CanvasXpress configuration passed to
#'   [canvasXpress::canvasXpress()].
#'
#' @return A `canvasXpress` htmlwidget.
#'
#' @examples
#' res <- data.frame(
#'     log2FoldChange = rnorm(100),
#'     padj           = runif(100),
#'     row.names      = paste0("gene", 1:100))
#' cxvolcano(res, title = "Volcano")
#'
#' @export
cxvolcano <- function(x, logfc = "log2FoldChange", pval = "padj", ...) {
    df <- as.data.frame(x)
    if (!all(c(logfc, pval) %in% colnames(df))) {
        stop(
            "Columns '", logfc, "' and '", pval,
            "' must both be present in 'x'.", call. = FALSE)
    }
    ok <- is.finite(df[[logfc]]) & is.finite(df[[pval]]) & df[[pval]] > 0
    df <- df[ok, , drop = FALSE]
    if (is.null(rownames(df))) {
        rownames(df) <- paste0("F", seq_len(nrow(df)))
    }

    mat <- t(as.matrix(data.frame(
        "log2FC"  = df[[logfc]],
        "negLogP" = -log10(df[[pval]]),
        row.names = rownames(df),
        check.names = FALSE)))

    canvasXpress::canvasXpress(
        data      = mat,
        graphType = "Scatter2D",
        ...)
}
