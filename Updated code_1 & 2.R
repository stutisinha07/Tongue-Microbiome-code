##############################################################################
# UPDATED CODE: keep the fuller dataset
# - Filter sparse genera only (>76% missing removed)
# - Keep ALL subjects with non-missing periogroup
# - Complete subject × genus grid
# - Fill absent combinations with 0
# - Add pseudocount +1
# - CLR transform
##############################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(haven)
  library(compositions)
  library(ggplot2)
  library(plotly)
  library(patchwork)
  library(mixOmics)
})

##############################################################################
# 0) HELPER FUNCTIONS
##############################################################################

# Count unique subjects per periogroup
count_perio <- function(df, id = "zz_nr", group = "periogroup", step = "step") {
  df %>%
    mutate(
      id2 = as.character(.data[[id]]),
      group2 = haven::as_factor(.data[[group]], levels = "labels")
    ) %>%
    distinct(id2, .keep_all = TRUE) %>%
    count(group2, name = "n_subjects", .drop = FALSE) %>%
    mutate(step = step, .before = 1) %>%
    rename(periogroup = group2)
}

# Build full subject x retained-genus matrix and CLR transform

build_full_clr_dataset <- function(df_long_raw,
                                   df_subjects,
                                   missing_cutoff = 76,
                                   id_col = "zz_nr",
                                   genus_col = "Genus",
                                   count_col = "counts") {
  
  dat <- df_long_raw %>%
    mutate(
      zz_nr = as.character(.data[[id_col]]),
      Genus = as.character(.data[[genus_col]]),
      counts = as.numeric(.data[[count_col]])
    )
  
  subs <- df_subjects %>%
    mutate(zz_nr = as.character(.data[[id_col]])) %>%
    distinct(zz_nr, .keep_all = TRUE)
  
  genus_missing_counts <- dat %>%
    mutate(counts_for_filter = na_if(counts, 0)) %>%
    group_by(Genus) %>%
    summarise(
      n_rows = n(),
      n_missing = sum(is.na(counts_for_filter)),
      pct_missing = round(100 * n_missing / n_rows, 2),
      .groups = "drop"
    ) %>%
    filter(pct_missing <= missing_cutoff) %>%
    arrange(desc(n_missing))
  
  kept_genera <- sort(unique(genus_missing_counts$Genus))
  
  dat_kept <- dat %>%
    filter(Genus %in% kept_genera) %>%
    mutate(
      Genus = stringr::str_replace_all(Genus, "\u00A0", " "),
      Genus = stringr::str_squish(Genus),
      Genus = if_else(is.na(Genus) | Genus == "", "UNKNOWN_GENUS", Genus)
    ) %>%
    group_by(zz_nr, Genus) %>%
    summarise(counts = sum(counts, na.rm = TRUE), .groups = "drop")
  
  kept_genera_clean <- sort(unique(dat_kept$Genus))
  
  X_wide <- tidyr::expand_grid(
    zz_nr = subs$zz_nr,
    Genus = kept_genera_clean
  ) %>%
    left_join(dat_kept, by = c("zz_nr", "Genus")) %>%
    mutate(counts = dplyr::coalesce(counts, 0)) %>%
    tidyr::pivot_wider(
      names_from   = Genus,
      values_from  = counts,
      names_prefix = "g_",
      names_repair = "unique"
    )
  
  mat <- as.matrix(dplyr::select(X_wide, -zz_nr)) + 1
  clr_mat <- compositions::clr(mat)
  rownames(clr_mat) <- X_wide$zz_nr
  
  X_clr_tbl <- dplyr::bind_cols(
    tibble::tibble(zz_nr = X_wide$zz_nr),
    tibble::as_tibble(clr_mat)
  )
  
  list(
    genus_missing_counts = genus_missing_counts,
    kept_genera = kept_genera_clean,
    X_wide = X_wide,
    Xmat_clr = clr_mat,
    X_clr_tbl = X_clr_tbl
  )
}

# Align periogroup to CLR matrix rows
make_pls_input <- function(Xmat_clr, df_meta) {
  meta_perio <- df_meta %>%
    dplyr::select(zz_nr, periogroup) %>%
    dplyr::mutate(
      zz_nr = as.character(zz_nr),
      periogroup = haven::as_factor(periogroup, levels = "labels")
    ) %>%
    dplyr::group_by(zz_nr) %>%
    dplyr::summarise(
      periogroup = dplyr::first(stats::na.omit(periogroup)),
      .groups = "drop"
    )
  
  ids_X <- rownames(Xmat_clr)
  
  align_df <- tibble::tibble(zz_nr = ids_X) %>%
    dplyr::left_join(meta_perio, by = "zz_nr")
  
  keep <- !is.na(align_df$periogroup)
  
  X_use <- Xmat_clr[keep, , drop = FALSE]
  y_use <- droplevels(factor(align_df$periogroup[keep]))
  
  list(
    meta_perio = meta_perio,
    align_df = align_df,
    keep = keep,
    X_use = X_use,
    y_use = y_use
  )
}

##############################################################################
# 1) SUBJECTS TO KEEP
#    Keep ALL subjects with non-missing periogroup
##############################################################################

subjects_s1 <- df_s1_bo_outcome %>%
  filter(!is.na(periogroup)) %>%
  distinct(zz_nr, .keep_all = TRUE)

subjects_s2 <- df_s2_bo_outcome %>%
  filter(!is.na(periogroup)) %>%
  distinct(zz_nr, .keep_all = TRUE)

cat("SHIP-1 subjects with non-missing periogroup:\n")
print(count_perio(subjects_s1, step = "subjects_s1"))

cat("\nSHIP-2 subjects with non-missing periogroup:\n")
print(count_perio(subjects_s2, step = "subjects_s2"))

##############################################################################
# 2) SHIP-1: sparse-genus filter + full subject grid + CLR
##############################################################################

ship1_res <- build_full_clr_dataset(
  df_long_raw = df_s1_bo_outcome_genus,
  df_subjects = subjects_s1,
  missing_cutoff = 76,
  id_col = "zz_nr",
  genus_col = "Genus",
  count_col = "counts"
)

genus_missing_counts_s1 <- ship1_res$genus_missing_counts
genera_kept_s1 <- ship1_res$kept_genera
X_s1 <- ship1_res$X_wide
Xmat_s1_clr <- ship1_res$Xmat_clr
X_s1_clr <- ship1_res$X_clr_tbl

cat("\nSHIP-1 number of retained genera:\n")
print(length(genera_kept_s1))

cat("\nSHIP-1 number of subjects in wide matrix:\n")
print(nrow(X_s1))

# Long CLR table for plots
df_s1_clr_long <- X_s1_clr %>%
  pivot_longer(
    cols = -zz_nr,
    names_to = "Genus",
    values_to = "clr_counts"
  )

##############################################################################
# 3) SHIP-2: sparse-genus filter + full subject grid + CLR
##############################################################################

ship2_res <- build_full_clr_dataset(
  df_long_raw = df_s2_bo_outcome_genus,
  df_subjects = subjects_s2,
  missing_cutoff = 76,
  id_col = "zz_nr",
  genus_col = "Genus",
  count_col = "counts"
)

genus_missing_counts_s2 <- ship2_res$genus_missing_counts
genera_kept_s2 <- ship2_res$kept_genera
X_s2 <- ship2_res$X_wide
Xmat_s2_clr <- ship2_res$Xmat_clr
X_s2_clr <- ship2_res$X_clr_tbl

cat("\nSHIP-2 number of retained genera:\n")
print(length(genera_kept_s2))

cat("\nSHIP-2 number of subjects in wide matrix:\n")
print(nrow(X_s2))

# Long CLR table for plots
df_s2_clr_long <- X_s2_clr %>%
  pivot_longer(
    cols = -zz_nr,
    names_to = "Genus",
    values_to = "clr_counts"
  )

##############################################################################
# 4) CHECK: did we keep the fuller dataset?
##############################################################################

ship1_final_subjects <- tibble(
  zz_nr = rownames(Xmat_s1_clr)
) %>%
  dplyr::left_join(
    subjects_s1 %>%
      dplyr::mutate(zz_nr = as.character(zz_nr)) %>%
      dplyr::select(zz_nr, periogroup),
    by = "zz_nr"
  )

ship2_final_subjects <- tibble(
  zz_nr = rownames(Xmat_s2_clr)
) %>%
  dplyr::left_join(
    subjects_s2 %>%
      dplyr::mutate(zz_nr = as.character(zz_nr)) %>%
      dplyr::select(zz_nr, periogroup),
    by = "zz_nr"
  )

count_perio(ship2_final_subjects, step = "SHIP2_full_grid_CLR")

cat("\nSHIP-1 counts after full-grid CLR build:\n")
print(count_perio(ship1_final_subjects, step = "SHIP1_full_grid_CLR"))

cat("\nSHIP-2 counts after full-grid CLR build:\n")
print(count_perio(ship2_final_subjects, step = "SHIP2_full_grid_CLR"))

subjects_s1 <- df_s1_bo_outcome %>%
  dplyr::filter(!is.na(periogroup)) %>%
  dplyr::mutate(zz_nr = as.character(zz_nr)) %>%
  dplyr::distinct(zz_nr, .keep_all = TRUE)

subjects_s2 <- df_s2_bo_outcome %>%
  dplyr::filter(!is.na(periogroup)) %>%
  dplyr::mutate(zz_nr = as.character(zz_nr)) %>%
  dplyr::distinct(zz_nr, .keep_all = TRUE)
##############################################################################
# 5) PLS-DA INPUTS
#    This should now keep all subjects with non-missing periogroup
##############################################################################

pls_input_s1 <- make_pls_input(Xmat_s1_clr, subjects_s1)
X_use <- pls_input_s1$X_use
y_use <- pls_input_s1$y_use
align_df <- pls_input_s1$align_df

pls_input_s2 <- make_pls_input(Xmat_s2_clr, subjects_s2)
X_use_s2 <- pls_input_s2$X_use
y_use_s2 <- pls_input_s2$y_use
align_df2 <- pls_input_s2$align_df

tibble::tibble(
  zz_nr = rownames(X_use),
  periogroup = y_use
) %>%
  dplyr::count(periogroup, name = "n_subjects")

tibble::tibble(
  zz_nr = rownames(X_use_s2),
  periogroup = y_use_s2
) %>%
  dplyr::count(periogroup, name = "n_subjects")


##############################################################################
# 6) OPTIONAL: PCA
##############################################################################

# -------------------
# SHIP-1 PCA
# -------------------
pca_res_s1 <- prcomp(Xmat_s1_clr, center = TRUE, scale. = FALSE)

ve_s1  <- pca_res_s1$sdev^2 / sum(pca_res_s1$sdev^2)
pct_s1 <- round(100 * ve_s1, 1)
pc_lab_s1 <- function(i) paste0("PC", i, " (", pct_s1[i], "%)")

pca_scores_3_s1 <- as.data.frame(pca_res_s1$x[, 1:3])
pca_scores_3_s1$Patient <- rownames(pca_scores_3_s1)

p1_s1 <- ggplot(pca_scores_3_s1, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  labs(
    title = "SHIP-1: PC1 vs PC2",
    subtitle = paste0("Combined variance: ", round(pct_s1[1] + pct_s1[2], 1), "%"),
    x = pc_lab_s1(1), y = pc_lab_s1(2)
  ) +
  theme_minimal()

p2_s1 <- ggplot(pca_scores_3_s1, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.7, color = "darkorange") +
  labs(
    title = "SHIP-1: PC1 vs PC3",
    subtitle = paste0("Combined variance: ", round(pct_s1[1] + pct_s1[3], 1), "%"),
    x = pc_lab_s1(1), y = pc_lab_s1(3)
  ) +
  theme_minimal()

p3_s1 <- ggplot(pca_scores_3_s1, aes(x = PC2, y = PC3)) +
  geom_point(alpha = 0.7, color = "seagreen") +
  labs(
    title = "SHIP-1: PC2 vs PC3",
    subtitle = paste0("Combined variance: ", round(pct_s1[2] + pct_s1[3], 1), "%"),
    x = pc_lab_s1(2), y = pc_lab_s1(3)
  ) +
  theme_minimal()

(p1_s1 | p2_s1) / p3_s1 +
  plot_annotation(
    title = "PCA on CLR-transformed genera (SHIP-1, fuller dataset)",
    caption = "All subjects with non-missing periogroup retained."
  )

# -------------------
# SHIP-2 PCA
# -------------------
pca_res_s2 <- prcomp(Xmat_s2_clr, center = TRUE, scale. = FALSE)

ve_s2  <- pca_res_s2$sdev^2 / sum(pca_res_s2$sdev^2)
pct_s2 <- round(100 * ve_s2, 1)
pc_lab_s2 <- function(i) paste0("PC", i, " (", pct_s2[i], "%)")

pca_scores_3_s2 <- as.data.frame(pca_res_s2$x[, 1:3])
pca_scores_3_s2$Patient <- rownames(pca_scores_3_s2)

p1_s2 <- ggplot(pca_scores_3_s2, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  labs(
    title = "SHIP-2: PC1 vs PC2",
    subtitle = paste0("Combined variance: ", round(pct_s2[1] + pct_s2[2], 1), "%"),
    x = pc_lab_s2(1), y = pc_lab_s2(2)
  ) +
  theme_minimal()

p2_s2 <- ggplot(pca_scores_3_s2, aes(x = PC1, y = PC3)) +
  geom_point(alpha = 0.7, color = "darkorange") +
  labs(
    title = "SHIP-2: PC1 vs PC3",
    subtitle = paste0("Combined variance: ", round(pct_s2[1] + pct_s2[3], 1), "%"),
    x = pc_lab_s2(1), y = pc_lab_s2(3)
  ) +
  theme_minimal()

p3_s2 <- ggplot(pca_scores_3_s2, aes(x = PC2, y = PC3)) +
  geom_point(alpha = 0.7, color = "seagreen") +
  labs(
    title = "SHIP-2: PC2 vs PC3",
    subtitle = paste0("Combined variance: ", round(pct_s2[2] + pct_s2[3], 1), "%"),
    x = pc_lab_s2(2), y = pc_lab_s2(3)
  ) +
  theme_minimal()

(p1_s2 | p2_s2) / p3_s2 +
  plot_annotation(
    title = "PCA on CLR-transformed genera (SHIP-2, fuller dataset)",
    caption = "All subjects with non-missing periogroup retained."
  )

##############################################################################
# 7) OPTIONAL: PLS-DA
##############################################################################

# -------------------
# SHIP-1 PLS-DA
# -------------------
ncomp_pls <- 3

plsda_s1 <- mixOmics::plsda(
  X = X_use,
  Y = y_use,
  ncomp = ncomp_pls
)

set.seed(123)
perf_plsda_s1 <- mixOmics::perf(
  plsda_s1,
  validation = "Mfold",
  folds = 5,
  progressBar = FALSE,
  nrepeat = 10
)

cat("\nSHIP-1 PLS-DA misclassification error:\n")
print(perf_plsda_s1$error.rate)

mixOmics::plotIndiv(
  plsda_s1,
  comp = c(1, 2),
  group = y_use,
  legend = TRUE,
  ellipse = TRUE,
  ind.names = FALSE,
  title = "PLS-DA (SHIP-1, fuller dataset): comp1 vs comp2"
)

mixOmics::plotIndiv(
  plsda_s1,
  comp = c(1, 3),
  group = y_use,
  legend = TRUE,
  ellipse = TRUE,
  ind.names = FALSE,
  title = "PLS-DA (SHIP-1, fuller dataset): comp1 vs comp3"
)

# -------------------
# SHIP-2 PLS-DA
# -------------------
plsda_s2 <- mixOmics::plsda(
  X = X_use_s2,
  Y = y_use_s2,
  ncomp = ncomp_pls
)

set.seed(123)
perf_plsda_s2 <- mixOmics::perf(
  plsda_s2,
  validation = "Mfold",
  folds = 5,
  progressBar = FALSE,
  nrepeat = 10
)

cat("\nSHIP-2 PLS-DA misclassification error:\n")
print(perf_plsda_s2$error.rate)

mixOmics::plotIndiv(
  plsda_s2,
  comp = c(1, 2),
  group = y_use_s2,
  legend = TRUE,
  ellipse = TRUE,
  ind.names = FALSE,
  title = "PLS-DA (SHIP-2, fuller dataset): comp1 vs comp2"
)

mixOmics::plotIndiv(
  plsda_s2,
  comp = c(1, 3),
  group = y_use_s2,
  legend = TRUE,
  ellipse = TRUE,
  ind.names = FALSE,
  title = "PLS-DA (SHIP-2, fuller dataset): comp1 vs comp3"
)

##############################################################################
# 8) FINAL CHECKS
##############################################################################

cat("\n====================================================\n")
cat("FINAL SUBJECT COUNTS USED IN PLS-DA\n")
cat("====================================================\n")

cat("\nSHIP-1:\n")
print(
  tibble(zz_nr = rownames(X_use), periogroup = y_use) %>%
    count(periogroup, name = "n_subjects")
)

cat("\nSHIP-2:\n")
print(
  tibble(zz_nr = rownames(X_use_s2), periogroup = y_use_s2) %>%
    count(periogroup, name = "n_subjects")
)

cat("\nSHIP-1 number of genera used:\n")
print(ncol(X_use))

cat("\nSHIP-2 number of genera used:\n")
print(ncol(X_use_s2))

