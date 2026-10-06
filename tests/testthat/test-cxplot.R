make_se <- function() {
    set.seed(1)
    SummarizedExperiment::SummarizedExperiment(
        assays  = list(counts = matrix(rnorm(40), nrow = 10,
                       dimnames = list(paste0("g", 1:10), paste0("s", 1:4)))),
        colData = S4Vectors::DataFrame(group = rep(c("A", "B"), 2)),
        rowData = S4Vectors::DataFrame(chr = sample(c("1", "2"), 10, TRUE)))
}

test_that("cxplot returns a canvasXpress widget for a SummarizedExperiment", {
    w <- cxplot(make_se())
    expect_s3_class(w, "canvasXpress")
    expect_s3_class(w, "htmlwidget")
})

test_that("cxplot carries assay dimensions and annotations", {
    w <- cxplot(make_se())
    expect_equal(length(w$x$data$y$vars), 10)
    expect_equal(length(w$x$data$y$smps), 4)
    # colData -> sample annotations are present
    expect_true(!is.null(w$x$data$x))
})

test_that("cxplot n= keeps only the most variable features", {
    w <- cxplot(make_se(), n = 5)
    expect_equal(length(w$x$data$y$vars), 5)
})

test_that("cxplot passes graphType through", {
    w <- cxplot(make_se(), graphType = "Dotplot")
    expect_equal(w$x$config$graphType, "Dotplot")
})

test_that("cxvolcano plots one point per gene", {
    res <- data.frame(log2FoldChange = rnorm(30), padj = runif(30),
                      row.names = paste0("gene", 1:30))
    v <- cxvolcano(res)
    expect_s3_class(v, "canvasXpress")
    # Genes are the variables (one point each); the fold-change and p-value are
    # the sample columns that drive the axes. (Transposing these collapses the
    # whole table to two points — the bug this guards against.)
    expect_setequal(unlist(v$x$data$y$vars), paste0("gene", 1:30))
    expect_setequal(unlist(v$x$data$y$smps), c("log2FC", "negLogP"))
})

test_that("cxvolcano errors on missing columns", {
    expect_error(cxvolcano(data.frame(a = 1, b = 2)), "must both be present")
})

test_that("cxvolcano labels only the requested genes", {
    res <- data.frame(log2FoldChange = rnorm(30), padj = runif(30),
                      row.names = paste0("gene", 1:30))
    v <- cxvolcano(res, label_genes = c("gene1", "gene2"))
    expect_equal(v$x$config$labelBy, "vars")
    expect_true("Highlight" %in% names(v$x$data$z))
    expect_equal(sum(unlist(v$x$data$z$Highlight) == "yes"), 2)
})

test_that("cxplot overlays and clustering reach the config", {
    w <- cxplot(make_se(), overlays = "group", var_overlays = "chr",
                cluster = TRUE)
    expect_equal(unlist(w$x$config$smpOverlays), "group")
    expect_equal(unlist(w$x$config$varOverlays), "chr")
    expect_true(isTRUE(w$x$config$samplesClustered))
    expect_true(isTRUE(w$x$config$variablesClustered))
})

test_that("cxsurvival builds a Kaplan-Meier widget", {
    cd <- data.frame(time = c(5, 12, 18, 24, 7, 15, 22, 30),
                     status = c(1, 1, 0, 1, 1, 0, 1, 0),
                     arm = rep(c("A", "B"), each = 4))
    w <- cxsurvival(cd, time = "time", status = "status", group = "arm")
    expect_s3_class(w, "canvasXpress")
    expect_equal(w$x$config$graphType, "KaplanMeier")
    expect_equal(w$x$config$colorBy, "arm")
    expect_equal(length(w$x$data$y$vars), 8)
})

test_that("cxsurvival errors on missing columns", {
    expect_error(cxsurvival(data.frame(a = 1), time = "time", status = "status"),
                 "missing")
})

test_that("cxboxplot builds a connected boxplot", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(expr = matrix(rnorm(8), nrow = 1,
                       dimnames = list("IL6", paste0("s", 1:8)))),
        colData = S4Vectors::DataFrame(
            Subject   = rep(paste0("P", 1:4), each = 2),
            Condition = rep(c("Baseline", "Week 12"), 4),
            row.names = paste0("s", 1:8)))
    w <- cxboxplot(se, feature = "IL6", group = "Condition",
                   connect = "Subject")
    expect_s3_class(w, "canvasXpress")
    expect_equal(w$x$config$graphType, "Boxplot")
    expect_equal(w$x$config$connectBy, "Subject")
    expect_equal(length(w$x$data$y$smps), 8)
})

test_that("cxboxplot errors on an unknown feature", {
    se <- SummarizedExperiment::SummarizedExperiment(
        assays  = list(expr = matrix(rnorm(4), nrow = 1,
                       dimnames = list("IL6", paste0("s", 1:4)))),
        colData = S4Vectors::DataFrame(g = rep(c("A", "B"), 2),
                                       row.names = paste0("s", 1:4)))
    expect_error(cxboxplot(se, feature = "NOPE", group = "g"),
                 "not a row of the assay")
})
