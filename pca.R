#!/usr/bin/env Rscript
# PCA module (BiocSingular-backed) for omnibenchmark.
#
# BiocSingular::runSVD dispatches on a BSPARAM object, so one code path covers
# several solvers and the solver becomes a benchmark axis rather than a module.
# Split out of the sibling `scrapper` module, which ships libscran's own PCA and
# has no business carrying a second SVD library.
#
# Implementation notes
# --------------------
# - runSVD expects cells-as-rows, so the gene-by-cell stage input is transposed
#   here; center = TRUE then centers per gene, the standard PCA convention.
# - Total variance is computed sparse-safe from row sums, not from the centered
#   matrix, so it is identical on every --dense / --deferred combination.

suppressPackageStartupMessages({
  library(Matrix)
  library(HDF5Array)
  library(BiocSingular)
  library(data.table)
})

# arg parsing
source("src/common/cli.R")
p <- arg_parser("PCA module (BiocSingular)")
p <- add_base_args(p)                    # --output_dir, --name
p <- add_stage_args(p, "PCA")     # the stage I/O contract
# your own method params — argparser directly (its add_argument requires `help`):
p <- add_argument(p, "--solver", type = "character", help = "exact, random or irlba")
p <- add_argument(p, "--n_components", type = "integer", help = "number of PCs")
p <- add_argument(p, "--random_seed", type = "integer", help = "seed")
# BLAS thread count. Snakemake exports OMP_NUM_THREADS = <rule threads>, which
# defaults to 1. omnibenchmark >= 7.0 should expose the number of threads if the module requires
# a number of cores.
# If you want to control it explicitely, RhpcBLASctl calls the
# library's own setter at runtime, which overrides the env var. 0 = leave alone.
p <- add_argument(p, "--blas_threads", type = "integer", default = 0L,
                  help = "BLAS threads (0 = inherit OMP_NUM_THREADS)")
# Dense vs sparse as its own axis: only the dense path puts real work through
# BLAS-3, and we want to compare the two directly.
p <- add_argument(p, "--dense", type = "character", default = "false",
                  help = "materialise the matrix dense before PCA (true/false)")
# Deferred centering. BiocSingular defaults every BSPARAM to deferred=FALSE,
# which materialises the centered matrix and destroys sparsity.
#
# Measured on BiocSingular 1.28.0, serial, fold=Inf (peak MB, FALSE/TRUE):
#   random 660/449   exact 652/642   irlba 358/359
# Only `random` gets a real effect. runExactSVD calls safe_svd(as.matrix(x)) so it
# materialises either way, and runIrlbaSVD's serial branch never reads `deferred` --
# it hands centering to irlba's own `center=` argument, which already defers.
# Both are kept in the grid on purpose: "deferral is not a per-method property" is
# the claim under test, and the null cells are what demonstrate it.
p <- add_argument(p, "--deferred", type = "character", default = "false",
                  help = "defer centering instead of materialising it (true/false)")
args <- parse_args(p)                    # argparser's own parser

# logging
cat(sprintf("Full command: %s\n", paste(commandArgs(trailingOnly = FALSE), collapse = " ")))
cat(sprintf("LOG: command line args\n----------------------------------\n"))
for (i in 1:length(args)) {
  cat(sprintf("  %s: %s\n", names(args)[i], args[[i]]))
}
cat(sprintf("----------------------------------\n"))


run_pca <- function(X, args) {
  # X: gene-by-cell sparse matrix (rows = genes).
  # Only --solver random consumes this: rsvd draws its test matrix from R's RNG.
  # exact is deterministic, and irlba's random start vector converges to `tol`
  # regardless, so the seed moves neither of them measurably.
  set.seed(args$random_seed)

  # Reference (netlib) BLAS is single-threaded by construction, so this is a
  # no-op there; this path is a true serial control.
  if (!is.na(args$blas_threads) && args$blas_threads > 0L) {
    RhpcBLASctl::blas_set_num_threads(args$blas_threads)
    cat(sprintf("LOG: blas threads set to %d (was %s)\n",
                args$blas_threads, Sys.getenv("OMP_NUM_THREADS", "unset")))
  }

  deferred <- identical(args$deferred, "true")
  bsparam <- switch(args$solver,
    random = RandomParam(deferred = deferred),
    exact  = ExactParam(deferred = deferred),
    irlba  = IrlbaParam(deferred = deferred),
    stop("unknown solver: ", args$solver)
  )
  cat(sprintf("LOG: centering %s\n",
              if (deferred) "deferred (operator)" else "materialised"))

  # runSVD expects cells-as-rows; center across genes (i.e. center=TRUE centers columns)
  svd <- runSVD(t(X), k = args$n_components, center = TRUE, BSPARAM = bsparam)
  embedding <- svd$u %*% diag(svd$d)        # (n_cells, n_components)
  loadings  <- svd$v                         # (n_genes, n_components)
  n         <- ncol(X)
  variance  <- svd$d^2 / (n - 1)
  # total variance: sum of per-gene variances, computed sparse-safe
  rs2 <- Matrix::rowSums(X^2)
  rs1 <- Matrix::rowSums(X)
  total_var <- sum((rs2 - rs1^2 / n) / (n - 1))

  variance_ratio <- variance / total_var
  # decorate embeddings/loadings w/ row/colnames
  rownames(embedding) <- colnames(X)
  colnames(embedding) <- paste0("PC", seq_len(ncol(embedding)))
  rownames(loadings)  <- rownames(X)
  colnames(loadings)  <- paste0("PC", seq_len(ncol(loadings)))

  # loadings, etc are here in case needed as output later
  list(
    embedding      = embedding,
    loadings       = loadings,
    variance       = as.double(variance),
    variance_ratio = as.double(variance_ratio)
  )
}


main <- function() {
  dir.create(args$output_dir, showWarnings = FALSE, recursive = TRUE)

  m <- TENxMatrix(args$normalized_selected_h5, group = "matrix")
  # Coerce straight from the DelayedArray to the target representation, so peak
  # memory reflects the format under test rather than a sparse copy plus a
  # dense one.
  if (identical(args$dense, "true")) {
    m <- as(m, "matrix")
    cat(sprintf("LOG: dense matrix: %.0f MB\n", as.numeric(object.size(m)) / 1e6))
  } else {
    m <- as(m, "dgCMatrix")
  }
  cat(sprintf("  matrix (genes x cells): %d x %d\n", nrow(m), ncol(m)))

  res <- run_pca(m, args)
  cat(sprintf("  embedding: %d x %d, loadings: %d x %d\n",
    nrow(res$embedding), ncol(res$embedding),
    nrow(res$loadings),  ncol(res$loadings)))

  out_embeddings_tsv <- file.path(args$output_dir, sprintf("%s_pcas.tsv", args$name))
  fwrite(data.frame(cell_id = rownames(res$embedding), res$embedding), out_embeddings_tsv,
    sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out_embeddings_tsv))

  out_loadings_tsv <- file.path(args$output_dir, sprintf("%s_loadings.tsv", args$name))
  fwrite(data.frame(gene = rownames(res$loadings), res$loadings), out_loadings_tsv,
    sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  wrote: %s\n", out_loadings_tsv))
}

if (sys.nframe() == 0L) {
  main()
}
