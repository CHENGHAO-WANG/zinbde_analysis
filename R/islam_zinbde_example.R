# Example ZINBDE analysis for the preprocessed Islam single-cell data.
# The script samples 100 expressed genes so it can demonstrate the workflow
# without running the full data set, then saves example outputs under output/.

required_packages <- c("zinbde", "future", "progressr", "glmmTMB")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    sprintf(
      "Install required package(s) before running this example: %s.",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

library(zinbde)
library(future)
library(progressr)

set.seed(20260426)

input_file <- file.path("data", "processed", "islam_dea_obj.rds")
output_dir <- file.path("output", "islam_zinbde_example")
wald_coef <- "cell_typestem_cell"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

islam <- readRDS(input_file)
cat("Original dimensions:", paste(dim(islam$count), collapse = " x "), "\n")

if (!"cell_type" %in% names(islam$metadata)) {
  stop(
    "Missing metadata column 'cell_type'. Run R/data_preprocessing.R before this example.",
    call. = FALSE
  )
}

keep_gene <- rowSums(islam$count > 0) >= 5L
islam_filtered <- islam[keep_gene, ]
cat("Filtered gene count:", nrow(islam_filtered$count), "\n")

islam_filtered <- estimate_size_factor(islam_filtered)

# randomly select 100 genes
sampled_genes <- sample(rownames(islam_filtered$count), size = 100L)
islam_example <- islam_filtered[sampled_genes, ]
cat("Sampled gene count:", nrow(islam_example$count), "\n")

# islam_example <- estimate_size_factor(islam_example)

workers <- max(1L, future::availableCores() - 1L)
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)

progressr::handlers("txtprogressbar")

fit_control_example <- fit_control(
  parallel = TRUE,
  progress = TRUE,
  future.chunk.size = 5L
)

glmmTMB_priors <- data.frame(
  prior = "cauchy(0,2.5)",
  class = "fixef_zi",
  component = "zi",
  # glmmTMB treats a blank fixed-effect coefficient as all non-intercept terms.
  # This avoids matching issues for formula-derived names that contain spaces.
  coef = c("(Intercept)", ""),
  stringsAsFactors = FALSE
)

cat("Fitting models with", workers, "future worker(s)...\n")
fit <- progressr::with_progress(
  fit_dea_model(
    object = islam_example,
    nb_formula = ~ cell_type,
    zi_formula = ~ cell_type,
    zero_inflation = "auto",
    nb_optimizer = "glmmTMB",
    zinb_optimizer = "glmmTMB",
    zi_prior = "none",
    control = fit_control_example,
    glmmTMB_args = list(
      priors = glmmTMB_priors,
      verbose = FALSE
    )
  )
)

cat("Running empirical Bayes dispersion shrinkage...\n")
fit_ebayes <- Ebayes(fit)

cat("Running Wald test for coefficient:", wald_coef, "\n")
wald <- wald_test(fit_ebayes, coef = wald_coef)

saveRDS(sampled_genes, file.path(output_dir, "selected_genes.rds"))
saveRDS(islam_example, file.path(output_dir, "islam_example_dea.rds"))
saveRDS(fit, file.path(output_dir, "islam_zinbde_fit.rds"))
saveRDS(fit_ebayes, file.path(output_dir, "islam_zinbde_ebayes.rds"))
saveRDS(wald, file.path(output_dir, "islam_wald_test.rds"))
utils::write.csv(
  wald,
  file = file.path(output_dir, "islam_wald_test.csv"),
  row.names = FALSE
)

cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat("Wald row count:", nrow(wald), "\n")
