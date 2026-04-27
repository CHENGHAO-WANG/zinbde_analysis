# Randomly split Islam fibroblast cells into equal-size groups, fit forced ZINB
# models on sampled genes, and compare Wald-test p-values before and after
# empirical Bayes dispersion shrinkage across repeated null splits.

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

set.seed(20260427)

input_file <- file.path("data", "processed", "islam_dea_obj.rds")
output_dir <- file.path("output", "islam_fibroblast_random_split_ebayes")
wald_results_file <- file.path(output_dir, "wald_p_values.rds")

tight_tol <- 1e-8
loose_tol <- 1e-2
em_max_iter <- 100L
component_max_iter <- 50L
sample_size <- as.integer(Sys.getenv("ISLAM_RANDOM_SPLIT_SAMPLE_SIZE", "100"))
n_iterations <- as.integer(Sys.getenv("ISLAM_RANDOM_SPLIT_ITERATIONS", "20"))
nb_formula <- ~ random_group
zi_formula <- ~ random_group
wald_coef <- "random_groupgroup_2"
plot_only <- tolower(Sys.getenv("ISLAM_RANDOM_SPLIT_PLOT_ONLY", "false")) %in%
  c("1", "true", "yes")

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

# Validate environment-variable overrides used for optional smoke tests.
as_positive_integer <- function(value, name) {
  if (length(value) != 1L || is.na(value) || value < 1L) {
    stop(sprintf("'%s' must be a positive integer.", name), call. = FALSE)
  }

  value
}

# Add iteration metadata to Wald-test rows and keep the output columns stable.
format_wald_rows <- function(wald, iteration, stage) {
  wald$iteration <- iteration
  wald$stage <- stage
  wald <- wald[
    ,
    c(
      "iteration",
      "gene",
      "stage",
      "model_type",
      "component",
      "term",
      "estimate",
      "se",
      "statistic",
      "p_value",
      "p_adj"
    )
  ]
  rownames(wald) <- NULL
  wald
}

# Save one histogram for one Wald-test metric, analysis stage, and component.
plot_histogram <- function(results, stage, component, value_column, filename, x_label) {
  plot_data <- results[
    results$stage == stage &
      results$component == component &
      is.finite(results[[value_column]]) &
      results[[value_column]] >= 0 &
      results[[value_column]] <= 1,
    ,
    drop = FALSE
  ]

  histogram <- ggplot(plot_data, aes(x = .data[[value_column]])) +
    geom_histogram(
      bins = 30,
      boundary = 0,
      closed = "left",
      color = "white",
      fill = "#2563EB"
    ) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
    labs(
      x = x_label,
      y = "Count",
      title = sprintf(
        "%s (%s, %s)",
        x_label,
        gsub("_", " ", stage, fixed = TRUE),
        component
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(
    filename = file.path(output_dir, filename),
    plot = histogram,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# Save p-value and adjusted p-value histograms for each Wald-test component.
plot_component_histograms <- function(results) {
  legacy_histograms <- file.path(
    output_dir,
    c(
      "p_value_before_ebayes_histogram.png",
      "p_value_after_ebayes_histogram.png",
      "p_adj_before_ebayes_histogram.png",
      "p_adj_after_ebayes_histogram.png"
    )
  )
  file.remove(legacy_histograms[file.exists(legacy_histograms)])

  expected_components <- c("nb", "zi", "joint")
  missing_components <- setdiff(expected_components, unique(results$component))

  if (length(missing_components) > 0L) {
    warning(
      sprintf(
        "Saved Wald results do not contain component(s): %s.",
        paste(missing_components, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  plot_specs <- expand.grid(
    stage = c("before_ebayes", "after_ebayes"),
    component = expected_components,
    value_column = c("p_value", "p_adj"),
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(plot_specs))) {
    spec <- plot_specs[i, ]
    plot_histogram(
      results = results,
      stage = spec$stage,
      component = spec$component,
      value_column = spec$value_column,
      filename = sprintf(
        "%s_%s_%s_histogram.png",
        spec$value_column,
        spec$stage,
        spec$component
      ),
      x_label = if (identical(spec$value_column, "p_value")) {
        "Wald p-value"
      } else {
        "Wald adjusted p-value"
      }
    )
  }
}

sample_size <- as_positive_integer(sample_size, "sample_size")
n_iterations <- as_positive_integer(n_iterations, "n_iterations")

if (plot_only) {
  if (!file.exists(wald_results_file)) {
    stop(
      sprintf("Cannot run in plot-only mode; missing %s.", wald_results_file),
      call. = FALSE
    )
  }

  wald_results <- readRDS(wald_results_file)
  plot_component_histograms(wald_results)
  cat("Loaded Wald results from:", normalizePath(wald_results_file, winslash = "/"), "\n")
  cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
  cat("Wald result rows:", nrow(wald_results), "\n")
  quit(save = "no")
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
cat("Initial filtered gene count:", nrow(islam_filtered$count), "\n")

islam_filtered <- estimate_size_factor(islam_filtered)
islam_fibroblast <- subset(islam_filtered, subset = cell_type == "fibroblast")
cat("Fibroblast cell count:", ncol(islam_fibroblast$count), "\n")

keep_fibroblast_gene <- rowSums(islam_fibroblast$count > 0) >= 5L &
  rowSums(islam_fibroblast$count == 0) >= 5L
islam_fibroblast <- islam_fibroblast[keep_fibroblast_gene, ]
cat("Fibroblast-filtered gene count:", nrow(islam_fibroblast$count), "\n")

if (nrow(islam_fibroblast$count) < sample_size) {
  stop(
    sprintf(
      "Only %d genes passed fibroblast filtering; %d are required for sampling.",
      nrow(islam_fibroblast$count),
      sample_size
    ),
    call. = FALSE
  )
}

cell_count <- ncol(islam_fibroblast$count)
if (cell_count < 2L || cell_count %% 2L != 0L) {
  stop(
    sprintf(
      "Expected an even number of fibroblast cells for equal splitting; found %d.",
      cell_count
    ),
    call. = FALSE
  )
}

sampled_genes <- sample(rownames(islam_fibroblast$count), size = sample_size)
islam_random_split <- islam_fibroblast[sampled_genes, ]
cat("Sampled gene count:", nrow(islam_random_split$count), "\n")

workers <- max(1L, future::availableCores() - 1L)
future::plan(future::multisession, workers = workers)
on.exit(future::plan(future::sequential), add = TRUE)

progressr::handlers("txtprogressbar")
options(progressr.enable = TRUE)

fit_control_random_split <- make_em_fisher_control()
fibroblast_cells <- colnames(islam_random_split$count)
half_cell_count <- length(fibroblast_cells) / 2L
wald_rows <- vector("list", n_iterations * 2L)
split_assignments <- vector("list", n_iterations)

cat(
  "Running",
  n_iterations,
  "random split iteration(s) with",
  workers,
  "future worker(s)...\n"
)

for (iteration in seq_len(n_iterations)) {
  cat("Iteration", iteration, "of", n_iterations, "\n")

  group_1_cells <- sample(fibroblast_cells, size = half_cell_count)
  random_group <- ifelse(fibroblast_cells %in% group_1_cells, "group_1", "group_2")

  iteration_object <- islam_random_split
  iteration_object$metadata$random_group <- factor(
    random_group,
    levels = c("group_1", "group_2")
  )
  rownames(iteration_object$metadata) <- fibroblast_cells

  split_assignments[[iteration]] <- data.frame(
    iteration = iteration,
    cell_id = fibroblast_cells,
    random_group = as.character(iteration_object$metadata$random_group),
    stringsAsFactors = FALSE
  )

  fit <- progressr::with_progress(
    fit_dea_model(
      object = iteration_object,
      nb_formula = nb_formula,
      zi_formula = zi_formula,
      zero_inflation = "zinb",
      nb_optimizer = "fisher",
      zinb_optimizer = "em_fisher",
      zi_prior = "jeffreys",
      control = fit_control_random_split
    )
  )

  wald_before <- wald_test(fit, coef = wald_coef)

  fit_ebayes <- Ebayes(fit)
  wald_after <- wald_test(fit_ebayes, coef = wald_coef)

  row_index <- (iteration - 1L) * 2L
  wald_rows[[row_index + 1L]] <- format_wald_rows(
    wald = wald_before,
    iteration = iteration,
    stage = "before_ebayes"
  )
  wald_rows[[row_index + 2L]] <- format_wald_rows(
    wald = wald_after,
    iteration = iteration,
    stage = "after_ebayes"
  )
}

wald_results <- do.call(rbind, wald_rows)
split_assignments <- do.call(rbind, split_assignments)

saveRDS(sampled_genes, file.path(output_dir, "selected_genes.rds"))
saveRDS(split_assignments, file.path(output_dir, "split_assignments.rds"))
saveRDS(wald_results, wald_results_file)
utils::write.csv(
  wald_results,
  file = file.path(output_dir, "wald_p_values.csv"),
  row.names = FALSE
)

plot_component_histograms(wald_results)

cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat("Wald result rows:", nrow(wald_results), "\n")
