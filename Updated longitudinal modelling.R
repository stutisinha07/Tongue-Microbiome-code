# ============================================================
# UPDATED LONGITUDINAL GENUS-WISE MODELLING
# Outcome: delta = SHIP-2 - SHIP-1 (CLR)
# Main predictor: baseline periogroup (Group 1 reference)
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
# 0) CHECK REQUIRED OBJECTS
# ------------------------------------------------------------
stopifnot(
  exists("delta_long"),
  exists("subjects_s1")
)

# ------------------------------------------------------------
# 1) HELPERS
# ------------------------------------------------------------
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

run_longitudinal_models_flexible <- function(dat, outcome_var = "delta") {
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

# ------------------------------------------------------------
# 2) BUILD BASELINE COVARIATE TABLE
# ------------------------------------------------------------
covars_s1_long <- build_covars_s1(subjects_s1)

cat("Baseline periogroup distribution:\n")
print(table(covars_s1_long$periogroup, useNA = "ifany"))

# ------------------------------------------------------------
# 3) MERGE DELTA OUTCOME WITH BASELINE COVARIATES
# ------------------------------------------------------------
mod_long <- delta_long %>%
  dplyr::mutate(
    zz_nr = as.character(zz_nr),
    Genus = as.character(Genus)
  ) %>%
  dplyr::left_join(covars_s1_long, by = "zz_nr") %>%
  dplyr::filter(!is.na(delta), !is.na(periogroup), !is.na(Genus))

cat("Number of genera to model:", dplyr::n_distinct(mod_long$Genus), "\n")
cat("Periogroup distribution in modelling data:\n")
print(table(mod_long$periogroup, useNA = "ifany"))

# ------------------------------------------------------------
# 4) FIT GENUS-WISE LONGITUDINAL MODELS
# ------------------------------------------------------------
cat("\n--- Running longitudinal models ---\n")
models_long <- run_longitudinal_models_flexible(mod_long, outcome_var = "delta")

if ("note" %in% names(models_long)) {
  cat("Skipped genera:\n")
  print(models_long %>% dplyr::filter(!is.na(note)))
}

# ------------------------------------------------------------
# 5) EXTRACT PERIOGROUP EFFECTS + FDR
# ------------------------------------------------------------
results_long <- extract_periogroup_results(models_long)

sig_long <- results_long %>% dplyr::filter(p_adj_fdr < 0.05)

cat("Significant genus x comparison hits (FDR < 0.05):", nrow(sig_long), "\n")
cat("Unique significant genera:", dplyr::n_distinct(sig_long$Genus), "\n")

print(sig_long)

write.csv(
  results_long,
  "SHIP_longitudinal_periogroup_genus_delta_models.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 6) OPTIONAL PLOTS
# ------------------------------------------------------------
if (nrow(results_long) > 0) {
  p_long <- ggplot(results_long, aes(x = estimate, y = -log10(p_adj_fdr))) +
    geom_point(alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    facet_wrap(~ comparison) +
    labs(
      title = "Baseline periogroup effects on genus delta CLR (SHIP-2 - SHIP-1)",
      x = "Effect size (delta CLR vs Group 1)",
      y = "-log10(FDR)"
    ) +
    theme_minimal(base_size = 12)
  
  print(p_long)
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
      x = "Effect size (delta CLR vs Group 1)",
      y = ""
    ) +
    theme_minimal(base_size = 12)
  
  print(p)
  invisible(top_df)
}

top10_long <- make_top10_forest(
  results_long,
  "Top 10 longitudinal periogroup effects on genus delta CLR"
)

save_results_table_pdf(
  results_long,
  "SHIP_longitudinal_periogroup_genus_results_paged.pdf",
  "Longitudinal periogroup genus regression results"
)

model_terms_long <- models_long %>%
  dplyr::filter(!is.na(Genus)) %>%
  dplyr::select(Genus, model_terms) %>%
  dplyr::distinct()

write.csv(
  model_terms_long,
  "SHIP_longitudinal_model_terms_by_genus.csv",
  row.names = FALSE
)