# Fit a forced ZINB model on sampled Islam genes, run empirical Bayes
# dispersion shrinkage, and plot original versus shrunk dispersion against
# mean normalized counts.

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
library(ggplot2)

set.seed(20260426)

input_file <- file.path("data", "processed", "islam_dea_obj.rds")
output_dir <- file.path("output", "islam_ebayes_dispersion_shrinkage")

tight_tol <- 1e-8
loose_tol <- 1e-2
em_max_iter <- 100L
component_max_iter <- 50L
sample_size <- 500L
nb_formula <- ~ cell_type
zi_formula <- ~ cell_type

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

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

cat("Input file:", normalizePath(input_file, winslash = "/", mustWork = TRUE), "\n")
islam <- readRDS(input_file)
cat("Original dimensions:", paste(dim(islam$count), collapse = " x "), "\n")

if (!"cell_type" %in% names(islam$metadata)) {
  stop(
    "Missing metadata column 'cell_type'. Run R/data_preprocessing.R before this script.",
    call. = FALSE
  )
}

keep_gene <- rowSums(islam$count > 0) >= 5L & rowSums(islam$count == 0) >= 5L
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
islam_ebayes_dispersion <- islam_filtered[sampled_genes, ]
cat("Sampled gene count:", nrow(islam_ebayes_dispersion$count), "\n")

workers <- max(1L, future::availableCores() - 1L)
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)

progressr::handlers("txtprogressbar")
options(progressr.enable = TRUE)

fit_control_ebayes <- make_em_fisher_control()

cat("Fitting forced ZINB models with", workers, "future worker(s)...\n")
fit <- progressr::with_progress(
  fit_dea_model(
    object = islam_ebayes_dispersion,
    nb_formula = nb_formula,
    zi_formula = zi_formula,
    zero_inflation = "zinb",
    nb_optimizer = "fisher",
    zinb_optimizer = "em_fisher",
    zi_prior = "jeffreys",
    control = fit_control_ebayes
  )
)

cat("Running empirical Bayes dispersion shrinkage...\n")
fit_ebayes <- Ebayes(fit)

if (
  is.null(fit_ebayes$ebayes$original_dispersion) ||
    is.null(fit_ebayes$ebayes$shrunk_dispersion)
) {
  stop("Ebayes output is missing original or shrunk dispersion values.", call. = FALSE)
}

gene_ids <- fit_ebayes$feature_ids
size_factor <- islam_ebayes_dispersion$metadata[colnames(islam_ebayes_dispersion$count), "size_factor"]
normalized_count <- sweep(
  islam_ebayes_dispersion$count[gene_ids, , drop = FALSE],
  MARGIN = 2L,
  STATS = size_factor,
  FUN = "/"
)

dispersion_shrinkage <- data.frame(
  gene = gene_ids,
  original_dispersion = unname(fit_ebayes$ebayes$original_dispersion[gene_ids]),
  shrunk_dispersion = unname(fit_ebayes$ebayes$shrunk_dispersion[gene_ids]),
  mean_normalized_count = unname(rowMeans(normalized_count)),
  stringsAsFactors = FALSE
)

plot_data <- rbind(
  data.frame(
    gene = dispersion_shrinkage$gene,
    mean_normalized_count = dispersion_shrinkage$mean_normalized_count,
    dispersion = dispersion_shrinkage$original_dispersion,
    dispersion_type = "Original",
    stringsAsFactors = FALSE
  ),
  data.frame(
    gene = dispersion_shrinkage$gene,
    mean_normalized_count = dispersion_shrinkage$mean_normalized_count,
    dispersion = dispersion_shrinkage$shrunk_dispersion,
    dispersion_type = "Shrunk",
    stringsAsFactors = FALSE
  )
)
plot_data <- plot_data[
  is.finite(plot_data$mean_normalized_count) &
    is.finite(plot_data$dispersion) &
    plot_data$mean_normalized_count > 0 &
    plot_data$dispersion > 0,
]
plot_data$dispersion_type <- factor(
  plot_data$dispersion_type,
  levels = c("Original", "Shrunk")
)
original_plot_data <- plot_data[plot_data$dispersion_type == "Original", ]

dispersion_plot <- ggplot(
  plot_data,
  aes(
    x = mean_normalized_count,
    y = dispersion,
    color = dispersion_type
  )
) +
  geom_point(alpha = 0.65, size = 1.5) +
  geom_smooth(
    data = original_plot_data,
    aes(x = mean_normalized_count, y = dispersion),
    method = "loess",
    formula = y ~ x,
    se = FALSE,
    color = "#111827",
    linewidth = 0.8,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = c("Original" = "#2563EB", "Shrunk" = "#DC2626"),
    name = "Dispersion"
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Mean normalized count (log10 scale)",
    y = "Dispersion (log10 scale)",
    title = "Original and shrunk dispersion by mean normalized count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

saveRDS(sampled_genes, file.path(output_dir, "selected_genes.rds"))
saveRDS(
  islam_ebayes_dispersion,
  file.path(output_dir, "islam_ebayes_dispersion_dea.rds")
)
saveRDS(fit, file.path(output_dir, "islam_zinb_fit.rds"))
saveRDS(fit_ebayes, file.path(output_dir, "islam_zinb_ebayes.rds"))
utils::write.csv(
  dispersion_shrinkage,
  file.path(output_dir, "dispersion_shrinkage.csv"),
  row.names = FALSE
)
ggsave(
  filename = file.path(output_dir, "dispersion_vs_mean_normalized_count.png"),
  plot = dispersion_plot,
  width = 7,
  height = 5,
  dpi = 300
)

cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat("Dispersion shrinkage rows:", nrow(dispersion_shrinkage), "\n")
