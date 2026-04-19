# ============================================================
# UPDATED sPLS-DA for the fuller CLR dataset
# Uses objects created by your updated preprocessing script:
#   subjects_s1, subjects_s2
#   Xmat_s1_clr, Xmat_s2_clr
#   make_pls_input()
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(haven)
  library(mixOmics)
  library(ggplot2)
})

# ------------------------------------------------------------
# 0) CHECK REQUIRED OBJECTS
# ------------------------------------------------------------
stopifnot(
  exists("Xmat_s1_clr"),
  exists("Xmat_s2_clr"),
  exists("subjects_s1"),
  exists("subjects_s2"),
  exists("make_pls_input")
)

# ------------------------------------------------------------
# 1) HELPER
# ------------------------------------------------------------
run_splsda_updated <- function(X, y, dataset_name = "SHIP-1",
                               ncomp = 3,
                               test.keepX = c(5, 10, 15, 20, 30, 45),
                               folds = 5,
                               nrepeat = 5) {
  
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"
  y <- droplevels(as.factor(y))
  
  max_ncomp <- max(1, min(ncomp, ncol(X), nlevels(y) - 1))
  if (max_ncomp < 1) max_ncomp <- 1
  
  keep_grid <- unique(sort(pmin(test.keepX, ncol(X))))
  keep_grid <- keep_grid[keep_grid >= 1]
  if (length(keep_grid) == 0) {
    keep_grid <- seq_len(min(5, ncol(X)))
  }
  
  cat("\n=========================\n")
  cat("sPLS-DA tuning —", dataset_name, "\n")
  cat("=========================\n")
  cat("n subjects:", nrow(X), "\n")
  cat("n genera:", ncol(X), "\n")
  cat("n classes:", nlevels(y), "\n")
  cat("components fitted:", max_ncomp, "\n")
  
  set.seed(123)
  tune_res <- mixOmics::tune.splsda(
    X = X,
    Y = y,
    ncomp = max_ncomp,
    test.keepX = keep_grid,
    validation = "Mfold",
    folds = folds,
    nrepeat = nrepeat,
    dist = "max.dist",
    measure = "BER",
    progressBar = FALSE
  )
  
  choice.keepX <- tune_res$choice.keepX
  choice.keepX <- as.numeric(choice.keepX)
  choice.keepX <- pmin(choice.keepX, ncol(X))
  
  if (length(choice.keepX) < max_ncomp) {
    choice.keepX <- c(
      choice.keepX,
      rep(tail(choice.keepX, 1), max_ncomp - length(choice.keepX))
    )
  }
  
  cat("\n[", dataset_name, "] Optimal keepX per component:\n", sep = "")
  print(choice.keepX)
  
  splsda_fit <- mixOmics::splsda(
    X = X,
    Y = y,
    ncomp = max_ncomp,
    keepX = choice.keepX,
    scale = FALSE
  )
  
  set.seed(123)
  perf_res <- mixOmics::perf(
    splsda_fit,
    validation = "Mfold",
    folds = folds,
    nrepeat = nrepeat,
    dist = "max.dist",
    auc = TRUE,
    progressBar = FALSE
  )
  
  cat("\n[", dataset_name, "] CV error rate:\n", sep = "")
  print(perf_res$error.rate)
  
  if (!is.null(perf_res$AUC)) {
    cat("\n[", dataset_name, "] AUC by component:\n", sep = "")
    print(perf_res$AUC)
  }
  
  mixOmics::plotIndiv(
    splsda_fit,
    comp = c(1, 2),
    group = y,
    legend = TRUE,
    ellipse = TRUE,
    ind.names = FALSE,
    title = paste0("sPLS-DA (", dataset_name, ", fuller dataset): comp1 vs comp2")
  )
  
  if (max_ncomp >= 3) {
    mixOmics::plotIndiv(
      splsda_fit,
      comp = c(1, 3),
      group = y,
      legend = TRUE,
      ellipse = TRUE,
      ind.names = FALSE,
      title = paste0("sPLS-DA (", dataset_name, ", fuller dataset): comp1 vs comp3")
    )
  }
  
  for (cc in seq_len(max_ncomp)) {
    mixOmics::plotLoadings(
      splsda_fit,
      comp = cc,
      method = "mean",
      contrib = "max",
      title = paste0("sPLS-DA loadings — ", dataset_name, " (comp ", cc, ")")
    )
  }
  
  invisible(list(
    fit = splsda_fit,
    tune = tune_res,
    perf = perf_res,
    keepX = choice.keepX
  ))
}

# ------------------------------------------------------------
# 2) ALIGN OUTCOMES TO CLR MATRICES
# ------------------------------------------------------------
pls_input_s1_spls <- make_pls_input(Xmat_s1_clr, subjects_s1)
X_s1_spls <- pls_input_s1_spls$X_use
y_s1_spls <- pls_input_s1_spls$y_use

pls_input_s2_spls <- make_pls_input(Xmat_s2_clr, subjects_s2)
X_s2_spls <- pls_input_s2_spls$X_use
y_s2_spls <- pls_input_s2_spls$y_use

cat("\nSHIP-1 final periogroup counts for sPLS-DA:\n")
print(table(y_s1_spls))

cat("\nSHIP-2 final periogroup counts for sPLS-DA:\n")
print(table(y_s2_spls))

# ------------------------------------------------------------
# 3) RUN sPLS-DA
# ------------------------------------------------------------
splsda_s1_res <- run_splsda_updated(
  X = X_s1_spls,
  y = y_s1_spls,
  dataset_name = "SHIP-1",
  ncomp = 3,
  test.keepX = c(5, 10, 15, 20, 30, 45),
  folds = 5,
  nrepeat = 5
)

splsda_fit_s1 <- splsda_s1_res$fit

splsda_s2_res <- run_splsda_updated(
  X = X_s2_spls,
  y = y_s2_spls,
  dataset_name = "SHIP-2",
  ncomp = 3,
  test.keepX = c(5, 10, 15, 20, 30, 45),
  folds = 5,
  nrepeat = 5
)

splsda_fit_s2 <- splsda_s2_res$fit

# Optional extras
plotLoadings(splsda_fit_s1, comp = 1)
plotVar(splsda_fit_s1, comp = c(1, 2))
plotVar(splsda_fit_s1, comp = c(1, 2), cutoff = 0.5)

plotLoadings(splsda_fit_s2, comp = 1)
plotVar(splsda_fit_s2, comp = c(1, 2))
plotVar(splsda_fit_s2, comp = c(1, 2), cutoff = 0.5)