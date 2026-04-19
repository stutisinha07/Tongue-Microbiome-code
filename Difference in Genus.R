# ============================================================
# UPDATED DIFFERENCE CODE
# Computes SHIP-2 - SHIP-1 on the CLR scale
# Works with either:
#   - all common genera
#   - common sPLS-DA-selected genera
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(mixOmics)
  library(haven)
})

# ------------------------------------------------------------
# 0) CHECK REQUIRED OBJECTS
# ------------------------------------------------------------
stopifnot(
  exists("Xmat_s1_clr"),
  exists("Xmat_s2_clr"),
  exists("subjects_s1")
)

# ------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------
standardize_genus <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("^g_", "", x)
  x
}

as_clr_matrix <- function(x, obj_name = deparse(substitute(x))) {
  if (is.data.frame(x) && "zz_nr" %in% names(x)) {
    rn <- as.character(x$zz_nr)
    x <- x %>% dplyr::select(-zz_nr)
    x <- as.data.frame(x, check.names = FALSE)
    rownames(x) <- rn
  }
  
  x <- as.data.frame(x, check.names = FALSE)
  
  if (is.null(rownames(x))) {
    stop(obj_name, " has no rownames. Rownames must be subject IDs.")
  }
  
  mat <- as.matrix(x)
  storage.mode(mat) <- "numeric"
  rownames(mat) <- as.character(rownames(mat))
  colnames(mat) <- standardize_genus(colnames(mat))
  mat
}

collapse_duplicate_cols <- function(mat) {
  cn <- colnames(mat)
  if (anyDuplicated(cn) == 0) return(mat)
  
  idx_list <- split(seq_along(cn), cn)
  
  out <- sapply(idx_list, function(idx) {
    if (length(idx) == 1) {
      mat[, idx]
    } else {
      rowMeans(mat[, idx, drop = FALSE], na.rm = TRUE)
    }
  })
  
  out <- as.matrix(out)
  rownames(out) <- rownames(mat)
  out
}

get_selected_genera <- function(fit, comps = 1:3) {
  comps <- comps[comps <= fit$ncomp]
  if (length(comps) == 0) return(character(0))
  
  out <- unique(unlist(lapply(comps, function(cc) {
    mixOmics::selectVar(fit, comp = cc)$name
  })))
  
  standardize_genus(out)
}

# ------------------------------------------------------------
# 2) PREPARE CLR MATRICES
# ------------------------------------------------------------
X1 <- as_clr_matrix(Xmat_s1_clr, "Xmat_s1_clr")
X2 <- as_clr_matrix(Xmat_s2_clr, "Xmat_s2_clr")

X1 <- collapse_duplicate_cols(X1)
X2 <- collapse_duplicate_cols(X2)

cat("SHIP-1 CLR matrix dimensions:", dim(X1), "\n")
cat("SHIP-2 CLR matrix dimensions:", dim(X2), "\n")

# ------------------------------------------------------------
# 3) CHOOSE GENERA SET
# ------------------------------------------------------------
use_spls_selected_genera <- FALSE

common_genera_all <- intersect(colnames(X1), colnames(X2))

if (use_spls_selected_genera) {
  stopifnot(exists("splsda_fit_s1"), exists("splsda_fit_s2"))
  
  sel1 <- get_selected_genera(splsda_fit_s1, comps = 1:3)
  sel2 <- get_selected_genera(splsda_fit_s2, comps = 1:3)
  
  selected_common <- intersect(sel1, sel2)
  
  if (length(selected_common) == 0) {
    stop("No common sPLS-DA-selected genera found between SHIP-1 and SHIP-2.")
  }
  
  genera_to_use <- selected_common
  cat("Using common sPLS-DA-selected genera:", length(genera_to_use), "\n")
} else {
  genera_to_use <- common_genera_all
  cat("Using all common genera:", length(genera_to_use), "\n")
}

# Optional overlap table
presence_tbl <- tibble(
  Genus = sort(unique(c(colnames(X1), colnames(X2))))
) %>%
  mutate(
    in_SHIP1 = Genus %in% colnames(X1),
    in_SHIP2 = Genus %in% colnames(X2),
    robustness = case_when(
      in_SHIP1 & in_SHIP2 ~ "Both",
      in_SHIP1 & !in_SHIP2 ~ "SHIP-1 only",
      !in_SHIP1 & in_SHIP2 ~ "SHIP-2 only",
      TRUE ~ NA_character_
    )
  )

print(presence_tbl %>% count(robustness))

# ------------------------------------------------------------
# 4) ALIGN SUBJECTS AND COMPUTE DELTA
# ------------------------------------------------------------
common_ids <- intersect(rownames(X1), rownames(X2))
common_ids <- sort(common_ids)

X1c <- X1[common_ids, genera_to_use, drop = FALSE]
X2c <- X2[common_ids, genera_to_use, drop = FALSE]

delta_clr <- X2c - X1c

delta_wide <- data.frame(
  zz_nr = common_ids,
  delta_clr,
  check.names = FALSE
)

delta_long <- delta_wide %>%
  as.data.frame() %>%
  tidyr::pivot_longer(
    cols = -zz_nr,
    names_to = "Genus",
    values_to = "delta"
  ) %>%
  dplyr::mutate(
    zz_nr = as.character(zz_nr),
    Genus = as.character(Genus)
  )

# ------------------------------------------------------------
# 5) SUBJECT-LEVEL SUMMARIES
# ------------------------------------------------------------
delta_summary <- delta_wide %>%
  dplyr::mutate(
    total_change = rowSums(dplyr::across(-zz_nr), na.rm = TRUE),
    abs_change   = rowSums(abs(as.matrix(dplyr::select(., -zz_nr))), na.rm = TRUE)
  ) %>%
  dplyr::select(zz_nr, total_change, abs_change)

head(delta_summary)

# ------------------------------------------------------------
# 6) GENUS-LEVEL SUMMARIES
# ------------------------------------------------------------
genus_change <- delta_long %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(
    mean_delta   = mean(delta, na.rm = TRUE),
    median_delta = median(delta, na.rm = TRUE),
    sd_delta     = sd(delta, na.rm = TRUE),
    n            = sum(!is.na(delta)),
    .groups      = "drop"
  )

cat("Number of genera in delta object:", nrow(genus_change), "\n")

ggplot(genus_change, aes(x = reorder(Genus, mean_delta), y = mean_delta)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Mean change in genus CLR values (SHIP-2 - SHIP-1)",
    x = "Genus",
    y = "Mean delta CLR"
  ) +
  theme_minimal(base_size = 12)

# ------------------------------------------------------------
# 7) DELTA BY BASELINE PERIOGROUP
# ------------------------------------------------------------
meta_perio_s1 <- subjects_s1 %>%
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

delta_long_meta <- delta_long %>%
  dplyr::left_join(meta_perio_s1, by = "zz_nr") %>%
  dplyr::filter(!is.na(periogroup)) %>%
  dplyr::mutate(periogroup = droplevels(as.factor(periogroup)))

ggplot(delta_long_meta, aes(x = periogroup, y = delta, fill = periogroup)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.alpha = 0.3) +
  facet_wrap(~ Genus, scales = "free_y") +
  labs(
    title = "Change in genus CLR values (SHIP-2 - SHIP-1) by baseline periogroup",
    x = "Baseline periogroup",
    y = "Delta CLR"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")
