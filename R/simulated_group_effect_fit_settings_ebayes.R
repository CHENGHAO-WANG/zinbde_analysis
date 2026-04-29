# Fit simulated group-effect data under selected ZINB settings, then compare
# Wald tests before and after mean-prior empirical Bayes dispersion shrinkage.

required_packages <- c("zinbde", "future", "progressr", "glmmTMB")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    sprintf(
      "Install required package(s) before running this script: %s.",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

library(zinbde)
library(future)
library(progressr)

set.seed(20260427)

input_file <- file.path(
  "output",
  "simulated_group_effect_count_data",
  "simulated_dea.rds"
)
output_dir <- file.path("output", "simulated_group_effect_fit_settings_ebayes")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tight_tol <- 1e-8
loose_tol <- 1e-2
em_max_iter <- 100L
component_max_iter <- 50L
normal_prior_sd <- c("(Intercept)" = 10, "group_2" = 2.5)
nb_formula <- ~ group_2
zi_formula <- ~ group_2
wald_coef <- "group_2"

# Build glmmTMB ZI priors with a separate intercept scale and shared scale for
# all non-intercept fixed effects.
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

# Build EM Fisher controls for the requested component and M-step tolerances,
# enabling Levenberg damping in the component fits.
make_em_fisher_control <- function(component_tol) {
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
        m_step_max_iter = 1L,
        m_step_tol = loose_tol
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

# Use the shared parallel/progress controls for glmmTMB fits.
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

# Add setting and stage metadata to Wald output while preserving test columns.
format_wald_rows <- function(wald, setting, stage) {
  wald$setting <- setting
  wald$stage <- stage
  wald <- wald[
    ,
    c(
      "setting",
      "stage",
      "gene",
      "model_type",
      "component",
      "term",
      "estimate",
      "se",
      "df",
      "statistic",
      "p_value",
      "p_adj"
    )
  ]
  rownames(wald) <- NULL
  wald
}

settings <- list(
  em_fisher__jeffreys__tight_components_loose_mstep = list(
    setting = "em_fisher__jeffreys__tight_components_loose_mstep",
    nb_optimizer = "fisher",
    zinb_optimizer = "em_fisher",
    zi_prior = "jeffreys",
    zi_prior_sd = NULL,
    control = make_em_fisher_control(component_tol = tight_tol),
    glmmTMB_args = NULL
  ),
  em_fisher__normal__tight_components_loose_mstep = list(
    setting = "em_fisher__normal__tight_components_loose_mstep",
    nb_optimizer = "fisher",
    zinb_optimizer = "em_fisher",
    zi_prior = "normal",
    zi_prior_sd = normal_prior_sd,
    control = make_em_fisher_control(component_tol = tight_tol),
    glmmTMB_args = NULL
  ),
  em_fisher__jeffreys__loose_components_loose_mstep = list(
    setting = "em_fisher__jeffreys__loose_components_loose_mstep",
    nb_optimizer = "fisher",
    zinb_optimizer = "em_fisher",
    zi_prior = "jeffreys",
    zi_prior_sd = NULL,
    control = make_em_fisher_control(component_tol = loose_tol),
    glmmTMB_args = NULL
  ),
  em_fisher__normal__loose_components_loose_mstep = list(
    setting = "em_fisher__normal__loose_components_loose_mstep",
    nb_optimizer = "fisher",
    zinb_optimizer = "em_fisher",
    zi_prior = "normal",
    zi_prior_sd = normal_prior_sd,
    control = make_em_fisher_control(component_tol = loose_tol),
    glmmTMB_args = NULL
  ),
  glmmTMB__cauchy = list(
    setting = "glmmTMB__cauchy",
    nb_optimizer = "glmmTMB",
    zinb_optimizer = "glmmTMB",
    zi_prior = "none",
    zi_prior_sd = NULL,
    control = make_glmmTMB_control(),
    glmmTMB_args = list(priors = make_glmmTMB_priors("cauchy"), verbose = FALSE)
  ),
  glmmTMB__normal = list(
    setting = "glmmTMB__normal",
    nb_optimizer = "glmmTMB",
    zinb_optimizer = "glmmTMB",
    zi_prior = "none",
    zi_prior_sd = NULL,
    control = make_glmmTMB_control(),
    glmmTMB_args = list(priors = make_glmmTMB_priors("normal"), verbose = FALSE)
  )
)

# Fit one setting with forced ZINB models for the simulated group coefficient.
run_setting <- function(object, setting) {
  args <- list(
    object = object,
    nb_formula = nb_formula,
    zi_formula = zi_formula,
    zero_inflation = "zinb",
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
simulated_input <- readRDS(input_file)
simulated_dea <- simulated_input[["object"]]

if (is.null(simulated_dea)) {
  stop("Input RDS must contain a DEA object at element 'object'.", call. = FALSE)
}
if (!"group_2" %in% names(simulated_dea$metadata)) {
  stop("Missing metadata column 'group_2' in simulated DEA object.", call. = FALSE)
}

cat("Original dimensions:", paste(dim(simulated_dea$count), collapse = " x "), "\n")
if ("size_factor" %in% names(simulated_dea$metadata)) {
  simulated_dea$metadata$size_factor <- NULL
  cat("Removed existing metadata column: size_factor\n")
}
simulated_dea <- estimate_size_factor(simulated_dea)
cat("Re-estimated size factors.\n")

workers <- max(1L, future::availableCores() - 1L)
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)

progressr::handlers("txtprogressbar")
options(progressr.enable = TRUE)

runtime_rows <- list()
fit_objects <- list()
fit_ebayes_objects <- list()
wald_rows <- vector("list", length(settings) * 2L)
wald_index <- 0L

cat(
  "Running",
  length(settings),
  "setting(s) with",
  workers,
  "future worker(s)...\n"
)

for (setting_name in names(settings)) {
  setting <- settings[[setting_name]]
  cat("Running setting:", setting_name, "\n")

  timing <- system.time({
    fit <- run_setting(simulated_dea, setting)
    wald_before <- wald_test(fit, coef = wald_coef)

    fit_ebayes <- Ebayes(fit, prior_type = "mean")
    wald_after <- wald_test(fit_ebayes, coef = wald_coef)
  })

  fit_objects[[setting_name]] <- fit
  fit_ebayes_objects[[setting_name]] <- fit_ebayes
  runtime_rows[[setting_name]] <- data.frame(
    setting = setting_name,
    elapsed_seconds = unname(timing[["elapsed"]]),
    user_seconds = unname(timing[["user.self"]]),
    system_seconds = unname(timing[["sys.self"]]),
    stringsAsFactors = FALSE
  )

  wald_index <- wald_index + 1L
  wald_rows[[wald_index]] <- format_wald_rows(
    wald = wald_before,
    setting = setting_name,
    stage = "before_ebayes"
  )
  wald_index <- wald_index + 1L
  wald_rows[[wald_index]] <- format_wald_rows(
    wald = wald_after,
    setting = setting_name,
    stage = "after_ebayes"
  )
}

runtime <- do.call(rbind, runtime_rows)
wald_results <- do.call(rbind, wald_rows)

saveRDS(simulated_dea, file.path(output_dir, "simulated_dea.rds"))
saveRDS(fit_objects, file.path(output_dir, "fit_objects.rds"))
saveRDS(fit_ebayes_objects, file.path(output_dir, "fit_ebayes_objects.rds"))
saveRDS(wald_results, file.path(output_dir, "wald_p_values.rds"))
utils::write.csv(
  wald_results,
  file = file.path(output_dir, "wald_p_values.csv"),
  row.names = FALSE
)
utils::write.csv(runtime, file.path(output_dir, "runtime.csv"), row.names = FALSE)

cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat("Completed settings:", length(settings), "\n")
cat("Wald result rows:", nrow(wald_results), "\n")
