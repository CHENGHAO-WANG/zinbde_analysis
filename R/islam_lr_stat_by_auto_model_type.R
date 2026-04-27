# Fit NB, ZINB, and auto zero-inflation models on sampled Islam genes.
# The script compares NB and ZINB log-likelihoods without running LR tests,
# then plots LR statistics grouped by the auto-selected model type.

required_packages <- c("zinbde", "future", "progressr", "ggplot2")
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

set.seed(20260426)

input_candidates <- file.path(
  "data",
  c("preprocessed", "processed"),
  "islam_dea_obj.rds"
)
input_file <- input_candidates[file.exists(input_candidates)][1L]
if (is.na(input_file)) {
  stop(
    sprintf(
      "Could not find the Islam DEA object. Checked: %s.",
      paste(input_candidates, collapse = ", ")
    ),
    call. = FALSE
  )
}

output_dir <- file.path("output", "islam_lr_stats_by_auto_model_type")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tight_tol <- 1e-8
loose_tol <- 1e-2
em_max_iter <- 100L
component_max_iter <- 50L
sample_size <- 200L
nb_formula <- ~ cell_type
zi_formula <- ~ cell_type

# Build the requested EM Fisher control: tight component convergence, loose
# M-step convergence, one M-step iteration, and Levenberg damping.
make_em_fisher_control <- function() {
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
        logistic_tol = tight_tol,
        use_levenberg = TRUE
      ),
      nb_control = nb_mean_control(
        nb_max_iter = component_max_iter,
        nb_tol = tight_tol,
        use_levenberg = TRUE
      ),
      dispersion_control = dispersion_control(
        cr_max_iter = component_max_iter,
        cr_tol = tight_tol
      )
    )
  )
}

# Run one zero-inflation mode with the shared optimizer and prior settings.
fit_zero_inflation_mode <- function(object, zero_inflation_mode, control) {
  progressr::with_progress(
    fit_dea_model(
      object = object,
      nb_formula = nb_formula,
      zi_formula = zi_formula,
      zero_inflation = zero_inflation_mode,
      nb_optimizer = "fisher",
      zinb_optimizer = "em_fisher",
      zi_prior = "jeffreys",
      control = control
    )
  )
}

# Convert a named logical vector from a fit object into gene-aligned values.
fit_flag <- function(fit, field, gene_ids) {
  values <- fit[[field]]
  if (is.null(values)) {
    return(rep(NA, length(gene_ids)))
  }
  unname(values[gene_ids])
}

cat("Input file:", normalizePath(input_file, winslash = "/", mustWork = TRUE), "\n")
islam <- readRDS(input_file)
cat("Original dimensions:", paste(dim(islam$count), collapse = " x "), "\n")

if (!"cell_type" %in% names(islam$metadata)) {
  stop(
    "Missing metadata column 'cell_type'. Run R/data_preprocessing.R before this script.",
    call. = FALSE
  )
}

keep_gene <- rowSums(islam$count > 0) >= 5L
islam_filtered <- islam[keep_gene, ]
cat("Filtered gene count:", nrow(islam_filtered$count), "\n")

if (nrow(islam_filtered$count) < sample_size) {
  stop(
    sprintf(
      "Only %d genes passed filtering; %d are required for sampling.",
      nrow(islam_filtered$count),
      sample_size
    ),
    call. = FALSE
  )
}

islam_filtered <- estimate_size_factor(islam_filtered)

sampled_genes <- sample(rownames(islam_filtered$count), size = sample_size)
islam_lr <- islam_filtered[sampled_genes, ]
cat("Sampled gene count:", nrow(islam_lr$count), "\n")

workers <- max(1L, future::availableCores() - 1L)
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)

progressr::handlers("txtprogressbar")
options(progressr.enable = TRUE)

fit_control_lr <- make_em_fisher_control()

cat("Fitting auto zero-inflation model with", workers, "future worker(s)...\n")
fit_auto <- fit_zero_inflation_mode(islam_lr, "auto", fit_control_lr)

cat("Fitting forced NB model...\n")
fit_nb <- fit_zero_inflation_mode(islam_lr, "nb", fit_control_lr)

cat("Fitting forced ZINB model...\n")
fit_zinb <- fit_zero_inflation_mode(islam_lr, "zinb", fit_control_lr)

gene_ids <- fit_auto$feature_ids
auto_model_types <- data.frame(
  gene_id = gene_ids,
  model_type = unname(fit_auto$model_type[gene_ids]),
  auto_converged = fit_flag(fit_auto, "converged", gene_ids),
  auto_fit_failed = fit_flag(fit_auto, "fit_failed", gene_ids),
  auto_optimizer = unname(fit_auto$optimizer[gene_ids]),
  stringsAsFactors = FALSE
)

lr_statistics <- data.frame(
  gene_id = gene_ids,
  auto_model_type = auto_model_types$model_type,
  logLik_nb = unname(fit_nb$logLik[gene_ids]),
  logLik_zinb = unname(fit_zinb$logLik[gene_ids]),
  nb_converged = fit_flag(fit_nb, "converged", gene_ids),
  nb_fit_failed = fit_flag(fit_nb, "fit_failed", gene_ids),
  zinb_converged = fit_flag(fit_zinb, "converged", gene_ids),
  zinb_fit_failed = fit_flag(fit_zinb, "fit_failed", gene_ids),
  auto_converged = auto_model_types$auto_converged,
  auto_fit_failed = auto_model_types$auto_fit_failed,
  stringsAsFactors = FALSE
)
lr_statistics$lr_stat <- 2 * (lr_statistics$logLik_zinb - lr_statistics$logLik_nb)

saveRDS(sampled_genes, file.path(output_dir, "selected_genes.rds"))
saveRDS(islam_lr, file.path(output_dir, "islam_lr_dea.rds"))
saveRDS(fit_auto, file.path(output_dir, "fit_auto.rds"))
saveRDS(fit_nb, file.path(output_dir, "fit_nb.rds"))
saveRDS(fit_zinb, file.path(output_dir, "fit_zinb.rds"))

utils::write.csv(
  auto_model_types,
  file.path(output_dir, "auto_model_types.csv"),
  row.names = FALSE
)
utils::write.csv(
  lr_statistics,
  file.path(output_dir, "lr_statistics.csv"),
  row.names = FALSE
)

plot_data <- lr_statistics[is.finite(lr_statistics$lr_stat), ]
plot_data$auto_model_type <- factor(plot_data$auto_model_type, levels = c("NB", "ZINB"))

lr_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = auto_model_type, y = lr_stat)
) +
  ggplot2::geom_boxplot(
    fill = "#93C5FD",
    color = "#1F2937",
    width = 0.55,
    outlier.shape = NA
  ) +
  ggplot2::geom_jitter(
    width = 0.15,
    height = 0,
    alpha = 0.65,
    size = 1.4,
    color = "#374151"
  ) +
  ggplot2::labs(
    x = "Auto-selected model type",
    y = "LR statistic: 2 * (logLik ZINB - logLik NB)",
    title = "LR statistics by auto-selected model type"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid.major.x = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = file.path(output_dir, "lr_statistics_by_auto_model_type.png"),
  plot = lr_plot,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat("LR statistic rows:", nrow(lr_statistics), "\n")
