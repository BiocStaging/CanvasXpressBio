# Bioconductor submission checklist — CanvasXpressBio

This package is a **bridge**: it depends on the CRAN `canvasXpress` package and adds
`SummarizedExperiment` visualization methods. It is a legitimate Bioconductor contribution
(visualization methods over Bioconductor core classes) and does **not** conflict with the CRAN
package, which stays on CRAN.

## Why a separate package (not submitting `canvasXpress` itself)

`canvasXpress` is already on CRAN and is general-purpose. Bioconductor does not accept a package
that is already on CRAN, and redirects general-purpose packages back to CRAN. A package must live
on CRAN **or** Bioconductor, not both. `CanvasXpressBio` is the Bioconductor-appropriate piece:
methods over `SummarizedExperiment` / `SingleCellExperiment` / `DESeqDataSet`.

## Pre-submission gate (run locally, must be clean)

```r
# From the package directory
devtools::document()                 # regenerate man/ + NAMESPACE
devtools::build_vignettes()
BiocManager::install("BiocCheck")    # not yet installed on the build box
BiocCheck::BiocCheck(".")            # must have 0 ERRORS
```

Also run `R CMD build .` then `R CMD check --no-manual <tarball>` (0 errors/warnings) and
`R CMD BiocCheck <tarball>`.

Known items to confirm before submission:
- [ ] `BiocCheck` returns 0 errors (install it first; it is not on the build box yet).
- [ ] Vignette builds with `BiocStyle` (listed by the vignette; add to `Suggests` if `BiocCheck`
      flags it — currently the vignette `output:` uses `BiocStyle::html_document`).
- [ ] Version stays `0.99.z` for the initial submission (Bioconductor requirement).
- [ ] All exported functions have runnable `@examples` (they do).
- [ ] Unit tests present and passing (`tests/testthat/`, currently green).

## Submission steps (maintainer — external, one-time)

1. Push `CanvasXpressBio/` to its **own public GitHub repository** (Bioconductor submits from a
   standalone repo, not a subdirectory). Tag/commit at the `0.99.0` version.
2. Open an issue at <https://github.com/Bioconductor/Contributions/issues> titled
   `CanvasXpressBio`, pasting the repo URL, per the New Package Submission template.
3. The Bioconductor Single Package Builder (SPB) runs `R CMD check` + `BiocCheck` on all
   platforms. Address any reviewer comments by pushing commits and bumping the `z` in `0.99.z`.
4. On acceptance the package is added to the next Bioconductor release; `git` access is granted
   for future updates (version continues on the Bioconductor even/odd devel/release scheme).

This step is outward-facing and is the maintainer's to perform; nothing here submits
automatically.
