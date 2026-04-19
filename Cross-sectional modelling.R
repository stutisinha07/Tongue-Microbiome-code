# ============================================================
# UPDATED CROSS-SECTIONAL genus-wise modelling at SHIP-1/SHIP-2
# Outcome: CLR-transformed genus abundance at one wave
# Main predictor: periogroup at the same wave (Group 1 reference)
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(haven)
  library(broom)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

# ------------------------------------------------------------
# 0) INPUT OBJECTS
# ------------------------------------------------------------
stopifnot(
  exists("Xmat_s1_clr"),
  exists("Xmat_s2_clr"),
  exists("subjects_s1"),
  exists("subjects_s2")
)

meta_s1_df <- subjects_s1
meta_s2_df <- subjects_s2

# TRUE  = use same genus set in both waves
# FALSE = use all genera separately in each wave
use_same_genera_both_waves <- TRUE

# If genera_to_use already exists from the difference script,
# it will be used when use_same_genera_both_waves = TRUE.

# ------------------------------------------------------------
# 1) HELPER FUNCTIONS
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
    x  <- x %>% dplyr::select(-zz_nr)
    x  <- as.data.frame(x, check.names = FALSE)
    rownames(x) <- rn
  }
  
  x <- as.data.frame(x, check.names = FALSE)
  
  if (is.null(rownames(x))) {
    stop(obj_name, " has no rownames. Rownames must be subject IDs (zz_nr).")
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

clr_to_long <- function(clr_mat, value_name = "clr_abundance") {
  as.data.frame(clr_mat, check.names = FALSE) %>%
    tibble::rownames_to_column("zz_nr") %>%
    tidyr::pivot_longer(
      cols = -zz_nr,
      names_to = "Genus",
      values_to = value_name
    ) %>%
    dplyr::mutate(
      zz_nr = as.character(zz_nr),
      Genus = as.character(Genus)
    )
}

extract_periogroup_results <- function(model_tbl) {
  model_tbl_clean <- model_tbl
  if ("note" %in% names(model_tbl_clean)) {
    model_tbl_clean <- model_tbl_clean %>% dplyr::filter(is.na(note))
  }
  
  model_tbl_clean %>%
    dplyr::filter(grepl("^periogroup", term)) %>%
    dplyr::mutate(
      p_adj_fdr = p.adjust(p.value, method = "fdr"),
      comparison = dplyr::case_when(
        term == "periogroup2" ~ "Group 2 vs Group 1",
        term == "periogroup3" ~ "Group 3 vs Group 1",
        term == "periogroup4" ~ "Group 4 vs Group 1",
        TRUE ~ term
      )
    ) %>%
    dplyr::select(
      Genus, comparison, estimate, std.error, statistic, p.value, p_adj_fdr, n
    ) %>%
    dplyr::arrange(p_adj_fdr)
}

run_cross_sectional_models_flexible <- function(dat, outcome_var = "clr_abundance") {
  genera <- sort(unique(dat$Genus))
  out <- vector("list", length(genera))
  names(out) <- genera
  
  candidate_terms <- c(
    "periogroup", "age", "sex", "smoking", "diabetes", "BMI", "hba1c",
    "dental_reason", "periodontal_tx", "power_brush", "interdental_use"
  )
  
  available_terms_global <- intersect(candidate_terms, names(dat))
  
  if (!"periogroup" %in% available_terms_global) {
    stop("periogroup is required in modelling data.")
  }
  
  for (g in genera) {
    df_g <- dat %>% dplyr::filter(Genus == g)
    
    needed_cols <- c(outcome_var, available_terms_global)
    df_g <- df_g %>%
      dplyr::filter(stats::complete.cases(df_g[, needed_cols, drop = FALSE])) %>%
      droplevels()
    
    if (nrow(df_g) < 20) {
      out[[g]] <- tibble::tibble(
        Genus = g,
        note = "skipped (too few complete-case rows)"
      )
      next
    }
    
    usable_terms <- c()
    
    for (v in available_terms_global) {
      x <- df_g[[v]]
      n_nonmiss_unique <- length(unique(x[!is.na(x)]))
      if (n_nonmiss_unique >= 2) {
        usable_terms <- c(usable_terms, v)
      }
    }
    
    if (!"periogroup" %in% usable_terms) {
      out[[g]] <- tibble::tibble(
        Genus = g,
        note = "skipped (periogroup has <2 levels after complete-case filtering)"
      )
      next
    }
    
    model_formula <- as.formula(
      paste(outcome_var, "~", paste(usable_terms, collapse = " + "))
    )
    
    fit <- lm(model_formula, data = df_g)
    
    out[[g]] <- broom::tidy(fit) %>%
      dplyr::mutate(
        Genus = g,
        n = nrow(df_g),
        model_terms = paste(usable_terms, collapse = " + ")
      )
  }
  
  dplyr::bind_rows(out)
}

build_covars_s1 <- function(df) {
  mit_cols   <- c("mit_3", "mit_4", "mit_5", "mit_6")
  whyza_cols <- paste0("whyza_", 1:10)
  
  required_cols <- c(
    "zz_nr", "periogroup", "sex", "age_ship1", "rau_1", "diab_1",
    "som_gew", "som_groe", "hba1c", "parobeh", "mit_2",
    mit_cols, whyza_cols
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("SHIP-1 metadata is missing these columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  df %>%
    dplyr::mutate(zz_nr = as.character(zz_nr)) %>%
    dplyr::distinct(zz_nr, .keep_all = TRUE) %>%
    dplyr::mutate(
      periogroup = haven::as_factor(periogroup, levels = "values"),
      periogroup = stats::relevel(periogroup, ref = "1"),
      sex = as.factor(sex),
      age = as.numeric(age_ship1),
      smoking = dplyr::case_when(
        rau_1 == 1 ~ "non-smoker",
        rau_1 == 2 ~ "former smoker",
        rau_1 == 3 ~ "current smoker",
        TRUE       ~ NA_character_
      ),
      smoking = factor(
        smoking,
        levels = c("non-smoker", "former smoker", "current smoker")
      ),
      diabetes = dplyr::case_when(
        diab_1 == 1 ~ "diabetic",
        diab_1 == 2 ~ "healthy",
        TRUE        ~ NA_character_
      ),
      diabetes = factor(diabetes, levels = c("healthy", "diabetic")),
      BMI = som_gew / ((som_groe / 100)^2),
      hba1c = as.numeric(hba1c),
      pain_reason = dplyr::if_any(dplyr::all_of(whyza_cols), ~ .x == 1),
      dental_reason = factor(
        dplyr::if_else(pain_reason, "pain/problem", "control"),
        levels = c("control", "pain/problem")
      ),
      periodontal_tx = as.factor(parobeh),
      power_brush = factor(
        dplyr::if_else(mit_2 == 1, "yes", "no"),
        levels = c("no", "yes")
      ),
      interdental_use = dplyr::if_any(dplyr::all_of(mit_cols), ~ .x == 1),
      interdental_use = factor(
        dplyr::if_else(interdental_use, "yes", "no"),
        levels = c("no", "yes")
      )
    ) %>%
    dplyr::select(
      zz_nr, periogroup, sex, age, smoking, diabetes, BMI, hba1c,
      dental_reason, periodontal_tx, power_brush, interdental_use
    )
}

build_covars_s2 <- function(df) {
  whyza_cols <- paste0("whyza", sprintf("%02d", 1:10))
  mit_cols   <- c("mit_03", "mit_04", "mit_05")
  
  required_cols <- c(
    "zz_nr", "periogroup", "SEX_SHIP2", "AGE_SHIP2",
    "som_gew", "som_groe", "parobeh", "mit_02",
    whyza_cols, mit_cols
  )
  
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("SHIP-2 metadata is missing these columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  df %>%
    dplyr::mutate(zz_nr = as.character(zz_nr)) %>%
    dplyr::distinct(zz_nr, .keep_all = TRUE) %>%
    dplyr::mutate(
      periogroup = haven::as_factor(periogroup, levels = "values"),
      periogroup = stats::relevel(periogroup, ref = "1"),
      sex = as.factor(SEX_SHIP2),
      age = as.numeric(AGE_SHIP2),
      BMI = som_gew / ((som_groe / 100)^2),
      pain_reason = dplyr::if_any(dplyr::all_of(whyza_cols), ~ .x == 1),
      dental_reason = factor(
        dplyr::if_else(pain_reason, "pain/problem", "control"),
        levels = c("control", "pain/problem")
      ),
      periodontal_tx = as.factor(parobeh),
      power_brush = factor(
        dplyr::if_else(mit_02 == 1, "yes", "no"),
        levels = c("no", "yes")
      ),
      interdental_use = dplyr::if_any(dplyr::all_of(mit_cols), ~ .x == 1),
      interdental_use = factor(
        dplyr::if_else(interdental_use, "yes", "no"),
        levels = c("no", "yes")
      )
    ) %>%
    dplyr::select(
      zz_nr, periogroup, sex, age, BMI,
      dental_reason, periodontal_tx, power_brush, interdental_use
    )
}

save_results_table_pdf <- function(results_df, file_name, title_text, rows_per_page = 25) {
  if (nrow(results_df) == 0) {
    message("No rows to save for ", file_name)
    return(invisible(NULL))
  }
  
  results_pdf <- results_df %>%
    dplyr::mutate(
      estimate   = round(estimate, 3),
      std.error  = round(std.error, 3),
      statistic  = round(statistic, 2),
      p.value    = signif(p.value, 3),
      p_adj_fdr  = signif(p_adj_fdr, 3)
    )
  
  n_pages <- ceiling(nrow(results_pdf) / rows_per_page)
  
  pdf(file_name, width = 11, height = 8.5)
  
  for (page in seq_len(n_pages)) {
    start <- (page - 1) * rows_per_page + 1
    end   <- min(page * rows_per_page, nrow(results_pdf))
    chunk <- results_pdf[start:end, ]
    
    grid.newpage()
    grid.text(
      paste0(title_text, " — page ", page, " / ", n_pages),
      x = 0.02, y = 0.98, just = c("left", "top"),
      gp = gpar(fontsize = 12, fontface = "bold")
    )
    
    gridExtra::grid.table(chunk, rows = NULL)
  }
  
  dev.off()
  invisible(NULL)
}

make_top10_forest <- function(results_df, plot_title) {
  if (nrow(results_df) == 0) {
    message("No results available for: ", plot_title)
    return(invisible(NULL))
  }
  
  top_df <- results_df %>%
    dplyr::mutate(
      conf.low  = estimate - 1.96 * std.error,
      conf.high = estimate + 1.96 * std.error
    ) %>%
    dplyr::arrange(p_adj_fdr) %>%
    dplyr::slice(1:min(10, n()))
  
  p <- ggplot(
    top_df,
    aes(
      x = estimate,
      y = reorder(paste(Genus, comparison, sep = " — "), estimate)
    )
  ) +
    geom_point() +
    geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(
      title = plot_title,
      x = "Effect size (CLR vs Group 1)",
      y = ""
    ) +
    theme_minimal(base_size = 12)
  
  print(p)
  invisible(top_df)
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
# 3) CHOOSE GENERA
# ------------------------------------------------------------
if (use_same_genera_both_waves) {
  if (exists("genera_to_use")) {
    genera_common <- intersect(intersect(colnames(X1), colnames(X2)), genera_to_use)
  } else {
    genera_common <- intersect(colnames(X1), colnames(X2))
  }
  
  if (length(genera_common) == 0) {
    stop("No common genera found between SHIP-1 and SHIP-2.")
  }
  
  X1_use <- X1[, genera_common, drop = FALSE]
  X2_use <- X2[, genera_common, drop = FALSE]
  
  cat("Using common genera in both waves:", length(genera_common), "\n")
} else {
  X1_use <- X1
  X2_use <- X2
  cat("Using all genera separately in each wave.\n")
}

# ------------------------------------------------------------
# 4) BUILD LONG OUTCOME TABLES
# ------------------------------------------------------------
s1_long <- clr_to_long(X1_use, value_name = "clr_abundance")
s2_long <- clr_to_long(X2_use, value_name = "clr_abundance")

# ------------------------------------------------------------
# 5) BUILD WAVE-SPECIFIC COVARIATE TABLES
# ------------------------------------------------------------
covars_s1 <- build_covars_s1(meta_s1_df)
covars_s2 <- build_covars_s2(meta_s2_df)

cat("SHIP-1 periogroup distribution:\n")
print(table(covars_s1$periogroup, useNA = "ifany"))

cat("SHIP-2 periogroup distribution:\n")
print(table(covars_s2$periogroup, useNA = "ifany"))

# ------------------------------------------------------------
# 6) MERGE OUTCOMES WITH COVARIATES
# ------------------------------------------------------------
mod_s1 <- s1_long %>%
  dplyr::left_join(covars_s1, by = "zz_nr") %>%
  dplyr::filter(!is.na(clr_abundance), !is.na(periogroup), !is.na(Genus))

mod_s2 <- s2_long %>%
  dplyr::left_join(covars_s2, by = "zz_nr") %>%
  dplyr::filter(!is.na(clr_abundance), !is.na(periogroup), !is.na(Genus))

cat("SHIP-1 genera to model:", dplyr::n_distinct(mod_s1$Genus), "\n")
cat("SHIP-2 genera to model:", dplyr::n_distinct(mod_s2$Genus), "\n")

# ------------------------------------------------------------
# 7) FIT GENUS-WISE CROSS-SECTIONAL MODELS
# ------------------------------------------------------------
cat("\n--- Running SHIP-1 cross-sectional models ---\n")
models_s1 <- run_cross_sectional_models_flexible(mod_s1, outcome_var = "clr_abundance")

cat("\n--- Running SHIP-2 cross-sectional models ---\n")
models_s2 <- run_cross_sectional_models_flexible(mod_s2, outcome_var = "clr_abundance")

if ("note" %in% names(models_s1)) {
  cat("Skipped SHIP-1 genera:\n")
  print(models_s1 %>% dplyr::filter(!is.na(note)))
}

if ("note" %in% names(models_s2)) {
  cat("Skipped SHIP-2 genera:\n")
  print(models_s2 %>% dplyr::filter(!is.na(note)))
}

# ------------------------------------------------------------
# 8) EXTRACT PERIOGROUP EFFECTS + FDR
# ------------------------------------------------------------
results_s1 <- extract_periogroup_results(models_s1)
results_s2 <- extract_periogroup_results(models_s2)

sig_s1 <- results_s1 %>% dplyr::filter(p_adj_fdr < 0.05)
sig_s2 <- results_s2 %>% dplyr::filter(p_adj_fdr < 0.05)

cat("SHIP-1 significant genus x comparison hits (FDR < 0.05):", nrow(sig_s1), "\n")
cat("SHIP-1 unique significant genera:", dplyr::n_distinct(sig_s1$Genus), "\n")

cat("SHIP-2 significant genus x comparison hits (FDR < 0.05):", nrow(sig_s2), "\n")
cat("SHIP-2 unique significant genera:", dplyr::n_distinct(sig_s2$Genus), "\n")

print(sig_s1)
print(sig_s2)

write.csv(results_s1, "SHIP1_cross_sectional_periogroup_genus_models.csv", row.names = FALSE)
write.csv(results_s2, "SHIP2_cross_sectional_periogroup_genus_models.csv", row.names = FALSE)

# ------------------------------------------------------------
# 9) OPTIONAL SUMMARY PLOTS
# ------------------------------------------------------------
if (nrow(results_s1) > 0) {
  p_s1 <- ggplot(results_s1, aes(x = estimate, y = -log10(p_adj_fdr))) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    facet_wrap(~ comparison) +
    labs(
      title = "SHIP-1 cross-sectional periogroup effects on genus CLR abundance",
      x = "Effect size (CLR vs Group 1)",
      y = "-log10(FDR)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p_s1)
}

if (nrow(results_s2) > 0) {
  p_s2 <- ggplot(results_s2, aes(x = estimate, y = -log10(p_adj_fdr))) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    facet_wrap(~ comparison) +
    labs(
      title = "SHIP-2 cross-sectional periogroup effects on genus CLR abundance",
      x = "Effect size (CLR vs Group 1)",
      y = "-log10(FDR)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p_s2)
}

top10_s1 <- make_top10_forest(
  results_s1,
  "Top 10 SHIP-1 periogroup effects on genus CLR abundance"
)

top10_s2 <- make_top10_forest(
  results_s2,
  "Top 10 SHIP-2 periogroup effects on genus CLR abundance"
)

save_results_table_pdf(
  results_s1,
  "SHIP1_cross_sectional_periogroup_genus_results_paged.pdf",
  "SHIP-1 cross-sectional periogroup genus regression results"
)

save_results_table_pdf(
  results_s2,
  "SHIP2_cross_sectional_periogroup_genus_results_paged.pdf",
  "SHIP-2 cross-sectional periogroup genus regression results"
)

model_terms_s1 <- models_s1 %>%
  dplyr::filter(!is.na(Genus)) %>%
  dplyr::select(Genus, model_terms) %>%
  dplyr::distinct()

model_terms_s2 <- models_s2 %>%
  dplyr::filter(!is.na(Genus)) %>%
  dplyr::select(Genus, model_terms) %>%
  dplyr::distinct()

write.csv(model_terms_s1, "SHIP1_cross_sectional_model_terms_by_genus.csv", row.names = FALSE)
write.csv(model_terms_s2, "SHIP2_cross_sectional_model_terms_by_genus.csv", row.names = FALSE)
