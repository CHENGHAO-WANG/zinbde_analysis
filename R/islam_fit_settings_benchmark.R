# Benchmark ZINBDE optimizer and zero-inflation-prior settings on the
# preprocessed Islam data. This script stops after fit_dea_model() and writes
# fit-time, zero-inflation decisions, coefficient, SE, and dispersion summaries.

required_packages <- c("zinbde", "future", "progressr", "glmmTMB")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    sprintf(
      "Install required package(s) before running this benchmark: %s.",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

library(zinbde)
library(future)
library(progressr)

set.seed(20260426)

input_candidates <- file.path(
  "data",
  c("preprocessed", "proprocessed", "processed"),
  "islam_dea_obj.rds"
)
input_file <- input_candidates[file.exists(input_candidates)][1L]
if (is.na(input_file)) {
  stop(
    sprintf(
      "Could not find the preprocessed Islam object. Checked: %s.",
      paste(input_candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

output_dir <- file.path("output", "islam_fit_settings_benchmark")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tight_tol <- 1e-8
loose_tol <- 1e-2
em_max_iter <- 100L
component_max_iter <- 50L
normal_prior_sd <- c("(Intercept)" = 10, "cell_typestem_cell" = 2.5)
nb_formula <- ~ cell_type
zi_formula <- ~ cell_type

make_glmmTMB_priors <- function(distribution) {
  if (!distribution %in% c("cauchy", "normal")) {
    stop("'distribution' must be 'cauchy' or 'normal'.", call. = FALSE)
  }

  data.frame(
    prior = c(
      sprintf("%s(0,10)", distribution),
      sprintf("%s(0,2.5)", distribution)
    ),
    class = "fixef_zi",
    component = "zi",
    # Blank coef means all non-intercept ZI fixed effects in glmmTMB priors.
    coef = c("(Intercept)", ""),
    stringsAsFactors = FALSE
  )
}

make_em_control <- function(component_tol, m_step_tol, m_step_max_iter) {
  fit_control(
    max_iter = em_max_iter,
    tol = tight_tol,
    parallel = TRUE,
    future.chunk.size = 5L,
    progress = TRUE,
    nb = nb_control(
      nb_mean_control = nb_mean_control(
        nb_max_iter = component_max_iter,
        nb_tol = tight_tol,
        use_levenberg = TRUE
      )
    ),
    zinb = zinb_control(
      em_fisher_control = zinb_em_fisher_control(
        max_iter = em_max_iter,
        tol = tight_tol
      ),
      m_step_control = zinb_mstep_control(
        m_step_max_iter = m_step_max_iter,
        m_step_tol = m_step_tol
      ),
      logistic_control = logistic_control(
        logistic_max_iter = component_max_iter,
        logistic_tol = component_tol,
        use_levenberg = TRUE
      ),
      nb_control = nb_mean_control(
        nb_max_iter = component_max_iter,
        nb_tol = component_tol,
        use_levenberg = TRUE
      ),
      dispersion_control = dispersion_control(
        cr_max_iter = component_max_iter,
        cr_tol = component_tol
      )
    )
  )
}

make_bfgs_control <- function() {
  fit_control(
    max_iter = em_max_iter,
    tol = tight_tol,
    parallel = TRUE,
    future.chunk.size = 5L,
    progress = TRUE,
    nb = nb_control(
      nb_mean_control = nb_mean_control(use_levenberg = TRUE),
      bfgs_control = nb_bfgs_control(reltol = tight_tol)
    ),
    zinb = zinb_control(
      bfgs_control = zinb_bfgs_control(reltol = tight_tol),
      logistic_control = logistic_control(use_levenberg = TRUE),
      nb_control = nb_mean_control(use_levenberg = TRUE)
    )
  )
}

make_glmmTMB_control <- function() {
  fit_control(
    max_iter = em_max_iter,
    tol = tight_tol,
    parallel = TRUE,
    future.chunk.size = 5L,
    progress = TRUE,
    nb = nb_control(
      nb_mean_control = nb_mean_control(use_levenberg = TRUE)
    ),
    zinb = zinb_control(
      logistic_control = logistic_control(use_levenberg = TRUE),
      nb_control = nb_mean_control(use_levenberg = TRUE)
    )
  )
}

em_variants <- list(
  em_tight_components_tight_mstep = list(
    component_tol = tight_tol,
    m_step_tol = tight_tol,
    m_step_max_iter = 50L
  ),
  gecm_tight_components_loose_mstep = list(
    component_tol = tight_tol,
    m_step_tol = loose_tol,
    m_step_max_iter = 1L
  ),
  ecm_loose_components_tight_mstep = list(
    component_tol = loose_tol,
    m_step_tol = tight_tol,
    m_step_max_iter = 50L
  ),
  gem_loose_components_loose_mstep = list(
    component_tol = loose_tol,
    m_step_tol = loose_tol,
    m_step_max_iter = 1L
  )
)

settings <- list()
for (prior in c("jeffreys", "normal")) {
  for (variant_name in names(em_variants)) {
    variant <- em_variants[[variant_name]]
    setting_name <- paste("em_fisher", prior, variant_name, sep = "__")
    settings[[setting_name]] <- list(
      setting = setting_name,
      nb_optimizer = "fisher",
      zinb_optimizer = "em_fisher",
      zi_prior = prior,
      zi_prior_sd = if (identical(prior, "normal")) normal_prior_sd else NULL,
      control = make_em_control(
        component_tol = variant$component_tol,
        m_step_tol = variant$m_step_tol,
        m_step_max_iter = variant$m_step_max_iter
      ),
      glmmTMB_args = NULL
    )
  }
}

settings$bfgs__jeffreys <- list(
  setting = "bfgs__jeffreys",
  nb_optimizer = "bfgs",
  zinb_optimizer = "bfgs",
  zi_prior = "jeffreys",
  zi_prior_sd = NULL,
  control = make_bfgs_control(),
  glmmTMB_args = NULL
)
settings$bfgs__normal <- list(
  setting = "bfgs__normal",
  nb_optimizer = "bfgs",
  zinb_optimizer = "bfgs",
  zi_prior = "normal",
  zi_prior_sd = normal_prior_sd,
  control = make_bfgs_control(),
  glmmTMB_args = NULL
)
settings$glmmTMB__cauchy <- list(
  setting = "glmmTMB__cauchy",
  nb_optimizer = "glmmTMB",
  zinb_optimizer = "glmmTMB",
  zi_prior = "none",
  zi_prior_sd = NULL,
  control = make_glmmTMB_control(),
  glmmTMB_args = list(priors = make_glmmTMB_priors("cauchy"), verbose = FALSE)
)
settings$glmmTMB__normal <- list(
  setting = "glmmTMB__normal",
  nb_optimizer = "glmmTMB",
  zinb_optimizer = "glmmTMB",
  zi_prior = "none",
  zi_prior_sd = NULL,
  control = make_glmmTMB_control(),
  glmmTMB_args = list(priors = make_glmmTMB_priors("normal"), verbose = FALSE)
)

coefficient_terms <- function(estimates, vcov_component) {
  coef_names <- names(estimates)
  if (length(coef_names) == length(estimates) && all(nzchar(coef_names))) {
    return(coef_names)
  }

  vcov_names <- if (is.matrix(vcov_component)) rownames(vcov_component) else NULL
  if (length(vcov_names) == length(estimates) && all(nzchar(vcov_names))) {
    return(vcov_names)
  }

  paste0("coef_", seq_along(estimates))
}

extract_standard_errors <- function(vcov_component, coef_names) {
  standard_errors <- rep(NA_real_, length(coef_names))
  names(standard_errors) <- coef_names

  if (!is.matrix(vcov_component)) {
    return(standard_errors)
  }

  variance <- diag(vcov_component)
  variance_names <- names(variance)
  if (length(variance_names) == length(variance) && all(nzchar(variance_names))) {
    matched <- intersect(variance_names, coef_names)
    standard_errors[matched] <- sqrt(pmax(variance[matched], 0))
    return(standard_errors)
  }

  if (length(variance) == length(coef_names)) {
    standard_errors[] <- sqrt(pmax(variance, 0))
  }

  standard_errors
}

coefficients_to_rows <- function(fit, setting) {
  rows <- list()

  for (gene_id in fit$feature_ids) {
    coef_entry <- fit$coefficients[[gene_id]]
    vcov_entry <- fit$vcov[[gene_id]]
    if (is.null(coef_entry)) {
      next
    }

    for (component in names(coef_entry)) {
      estimates <- coef_entry[[component]]
      if (is.null(estimates) || length(estimates) == 0L) {
        next
      }

      vcov_component <- if (is.null(vcov_entry)) NULL else vcov_entry[[component]]
      terms <- coefficient_terms(estimates, vcov_component)
      names(estimates) <- terms
      standard_errors <- extract_standard_errors(
        vcov_component = vcov_component,
        coef_names = terms
      )

      component_rows <- data.frame(
        setting = setting,
        gene_id = gene_id,
        component = component,
        term = terms,
        estimate = unname(estimates),
        std_error = unname(standard_errors[terms]),
        stringsAsFactors = FALSE
      )
      rows[[length(rows) + 1L]] <- component_rows
    }
  }

  if (length(rows) == 0L) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

zi_decisions_to_rows <- function(fit, setting) {
  do.call(
    rbind,
    lapply(fit$feature_ids, function(gene_id) {
      zi_diag <- fit$diagnostics[[gene_id]]$zero_inflation
      if (is.null(zi_diag)) {
        zi_diag <- list()
      }

      data.frame(
        setting = setting,
        gene_id = gene_id,
        model_type = unname(fit$model_type[gene_id]),
        needs_zero_inflation = identical(unname(fit$model_type[gene_id]), "ZINB"),
        decision_reason = if (is.null(zi_diag$decision_reason)) NA_character_ else zi_diag$decision_reason,
        test_performed = if (is.null(zi_diag$test_performed)) NA else zi_diag$test_performed,
        p_value = if (is.null(zi_diag$p_value)) NA_real_ else zi_diag$p_value,
        observed_zero_proportion = if (is.null(zi_diag$observed_zero_proportion)) NA_real_ else zi_diag$observed_zero_proportion,
        expected_zero_proportion = if (is.null(zi_diag$expected_zero_proportion)) NA_real_ else zi_diag$expected_zero_proportion,
        converged = unname(fit$converged[gene_id]),
        fit_failed = unname(fit$fit_failed[gene_id]),
        optimizer = unname(fit$optimizer[gene_id]),
        stringsAsFactors = FALSE
      )
    })
  )
}

dispersion_to_rows <- function(fit, setting) {
  data.frame(
    setting = setting,
    gene_id = names(fit$dispersion),
    dispersion = unname(fit$dispersion),
    stringsAsFactors = FALSE
  )
}

run_setting <- function(object, setting) {
  args <- list(
    object = object,
    nb_formula = nb_formula,
    zi_formula = zi_formula,
    zero_inflation = "auto",
    nb_optimizer = setting$nb_optimizer,
    zinb_optimizer = setting$zinb_optimizer,
    zi_prior = setting$zi_prior,
    control = setting$control,
    glmmTMB_args = setting$glmmTMB_args
  )

  if (!is.null(setting$zi_prior_sd)) {
    args$zi_prior_sd <- setting$zi_prior_sd
  }

  progressr::with_progress(do.call(fit_dea_model, args))
}

cat("Input file:", normalizePath(input_file, winslash = "/", mustWork = TRUE), "\n")
islam <- readRDS(input_file)
cat("Original dimensions:", paste(dim(islam$count), collapse = " x "), "\n")

if (!"cell_type" %in% names(islam$metadata)) {
  stop(
    "Missing metadata column 'cell_type'. Run R/data_preprocessing.R before this benchmark.",
    call. = FALSE
  )
}

keep_gene <- rowSums(islam$count > 0) >= 5L
islam_filtered <- islam[keep_gene, ]
cat("Filtered gene count:", nrow(islam_filtered$count), "\n")

islam_filtered <- estimate_size_factor(islam_filtered)

sampled_genes <- sample(rownames(islam_filtered$count), size = 100L)
islam_benchmark <- islam_filtered[sampled_genes, ]
cat("Sampled gene count:", nrow(islam_benchmark$count), "\n")

workers <- max(1L, future::availableCores() - 1L)
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)

progressr::handlers("txtprogressbar")
options(progressr.enable = TRUE)

runtime_rows <- list()
zi_decision_rows <- list()
coefficient_rows <- list()
dispersion_rows <- list()
fit_objects <- list()

for (setting_name in names(settings)) {
  setting <- settings[[setting_name]]
  cat("Running setting:", setting_name, "\n")

  timing <- system.time({
    fit <- run_setting(islam_benchmark, setting)
  })

  fit_objects[[setting_name]] <- fit
  runtime_rows[[setting_name]] <- data.frame(
    setting = setting_name,
    elapsed_seconds = unname(timing[["elapsed"]]),
    user_seconds = unname(timing[["user.self"]]),
    system_seconds = unname(timing[["sys.self"]]),
    stringsAsFactors = FALSE
  )
  zi_decision_rows[[setting_name]] <- zi_decisions_to_rows(fit, setting_name)
  coefficient_rows[[setting_name]] <- coefficients_to_rows(fit, setting_name)
  dispersion_rows[[setting_name]] <- dispersion_to_rows(fit, setting_name)
}

runtime <- do.call(rbind, runtime_rows)
zi_decisions <- do.call(rbind, zi_decision_rows)
coefficients <- do.call(rbind, coefficient_rows)
dispersion <- do.call(rbind, dispersion_rows)

saveRDS(sampled_genes, file.path(output_dir, "selected_genes.rds"))
saveRDS(islam_benchmark, file.path(output_dir, "islam_benchmark_dea.rds"))
saveRDS(fit_objects, file.path(output_dir, "fit_objects.rds"))

utils::write.csv(runtime, file.path(output_dir, "runtime.csv"), row.names = FALSE)
utils::write.csv(zi_decisions, file.path(output_dir, "zi_decisions.csv"), row.names = FALSE)
utils::write.csv(coefficients, file.path(output_dir, "coefficients.csv"), row.names = FALSE)
utils::write.csv(dispersion, file.path(output_dir, "dispersion.csv"), row.names = FALSE)

cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat("Completed settings:", length(settings), "\n")
