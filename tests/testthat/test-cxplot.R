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

test_that("cxvolcano builds a two-variable scatter", {
    res <- data.frame(log2FoldChange = rnorm(30), padj = runif(30),
                      row.names = paste0("gene", 1:30))
    v <- cxvolcano(res)
    expect_s3_class(v, "canvasXpress")
    expect_setequal(unlist(v$x$data$y$vars), c("log2FC", "negLogP"))
})

test_that("cxvolcano errors on missing columns", {
    expect_error(cxvolcano(data.frame(a = 1, b = 2)), "must both be present")
})
