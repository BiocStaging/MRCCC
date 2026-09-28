# SummarizedExperiment input support for MRCCC.

#' Extract one gene from a SummarizedExperiment as a donor-level vector
#'
#' Extracts a single gene (row) from a `SummarizedExperiment` and returns it
#' as a plain named numeric vector suitable for the `X`, `Z` or `Y` argument
#' of [mr_ccc()]. The expected layout is donor-level (pseudo-bulk)
#' expression with genes as rows and donors (samples) as columns, so the
#' returned vector holds one value per donor, named by the column names of
#' the object.
#'
#' @param se A `SummarizedExperiment` with genes as rows and donors
#'   (samples) as columns.
#' @param gene A single character string: the name of the gene to extract.
#'   Must match one of `rownames(se)`.
#' @param assay The assay to read from: either a single positive integer
#'   index or the name of an assay of `se`. Defaults to the first assay.
#'
#' @return A named numeric vector of length `ncol(se)` with names
#'   `colnames(se)`, holding the expression of `gene` in the selected assay.
#'
#' @seealso [mr_ccc()]
#'
#' @examples
#' sim <- simulate_mrccc(n = 60, seed = 1)
#' m <- rbind(LIG1 = as.numeric(sim$X), REC1 = as.numeric(sim$Z))
#' colnames(m) <- paste0("donor", seq_len(ncol(m)))
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(logcounts = m)
#' )
#' x <- se_row(se, "LIG1")
#' head(x)
#'
#' @export
se_row <- function(se, gene, assay = 1L) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  if (!is.character(gene) || length(gene) != 1L || is.na(gene)) {
    stop("'gene' must be a single character string.", call. = FALSE)
  }
  rn <- rownames(se)
  if (is.null(rn) || !(gene %in% rn)) {
    stop("Gene '", gene, "' was not found in rownames(se).", call. = FALSE)
  }
  a_names <- SummarizedExperiment::assayNames(se)
  n_assay <- length(SummarizedExperiment::assays(se))
  if (is.character(assay)) {
    if (length(assay) != 1L || is.na(assay) || !(assay %in% a_names)) {
      stop("'assay' must be the name of an assay of 'se'.", call. = FALSE)
    }
  } else {
    if (!is.numeric(assay) || length(assay) != 1L || !is.finite(assay) ||
        assay != round(assay) || assay < 1 || assay > n_assay) {
      stop("'assay' must be an assay name or a single integer between 1 ",
           "and ", n_assay, ".", call. = FALSE)
    }
    assay <- as.integer(assay)
  }
  row <- SummarizedExperiment::assay(se, assay)[gene, ]
  if (!is.numeric(row)) {
    stop("The row extracted from assay(se, assay) is not numeric.",
         call. = FALSE)
  }
  out <- as.numeric(row)
  names(out) <- colnames(se)
  out
}

#' Reduce a single-row SummarizedExperiment to a named numeric vector
#'
#' Applied to the `X`, `Z` and `Y` arguments of [mr_ccc()] before
#' validation. Objects that are not a `SummarizedExperiment` pass through
#' unchanged. A `SummarizedExperiment` with exactly one row is reduced to
#' its first assay row through [se_row()]; any other `SummarizedExperiment`
#' is an error that points the user to [se_row()].
#'
#' @param x Object supplied by the user.
#' @param arg Argument name used in error messages.
#' @return `x` unchanged, or a named numeric vector.
#' @noRd
#' @keywords internal
coerce_se_vector <- function(x, arg) {
  if (!methods::is(x, "SummarizedExperiment")) {
    return(x)
  }
  if (nrow(x) != 1L) {
    stop("'", arg, "' is a SummarizedExperiment with ", nrow(x), " rows; ",
         "use se_row() to select a single gene, for example ",
         arg, " = se_row(se, \"GENE\").", call. = FALSE)
  }
  if (is.null(rownames(x))) {
    rownames(x) <- "feature1"
  }
  se_row(x, rownames(x)[1L], assay = 1L)
}
