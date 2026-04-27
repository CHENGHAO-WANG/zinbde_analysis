# Simulate a compact ZINB count data set with known group effects.
# The script uses constant coefficients, size factors, and dispersions so the
# resulting truth table has four balanced gene categories.

required_packages <- c("zinbde")
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

seed <- 20260427
n_cells <- 50L
n_genes <- 100L
n_cells_per_group <- 25L
n_genes_per_category <- 25L

size_factor_value <- 1
dispersion_value <- 0.2
nb_intercept <- log(5)
nb_group_effect <- log(2)
zi_intercept <- qlogis(0.2)
zi_group_effect <- -1

output_dir <- file.path("output", "simulated_group_effect_count_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cell_id <- sprintf("cell_%03d", seq_len(n_cells))
gene_id <- sprintf("gene_%03d", seq_len(n_genes))
group <- rep(c("group_1", "group_2"), each = n_cells_per_group)
group_2 <- as.numeric(group == "group_2")

nb_covariates <- cbind(intercept = 1, group_2 = group_2)
zi_covariates <- cbind(intercept = 1, group_2 = group_2)
rownames(nb_covariates) <- cell_id
rownames(zi_covariates) <- cell_id

gene_category <- rep(
  c("null", "nb_only", "zi_only", "zi_and_nb"),
  each = n_genes_per_category
)

nb_parameters <- cbind(
  intercept = rep(nb_intercept, n_genes),
  group_2 = ifelse(gene_category %in% c("nb_only", "zi_and_nb"), nb_group_effect, 0)
)
zi_parameters <- cbind(
  intercept = rep(zi_intercept, n_genes),
  group_2 = ifelse(gene_category %in% c("zi_only", "zi_and_nb"), zi_group_effect, 0)
)
rownames(nb_parameters) <- gene_id
rownames(zi_parameters) <- gene_id

size_factor <- rep(size_factor_value, n_cells)
dispersion <- rep(dispersion_value, n_genes)
names(size_factor) <- cell_id
names(dispersion) <- gene_id

simulated <- simulate_dea_data(
  nb_covariates = nb_covariates,
  nb_parameters = nb_parameters,
  size_factor = size_factor,
  dispersion = dispersion,
  zi_covariates = zi_covariates,
  zi_parameters = zi_parameters,
  seed = seed
)

cell_metadata <- data.frame(
  cell_id = cell_id,
  group = group,
  group_2 = group_2,
  size_factor = size_factor,
  stringsAsFactors = FALSE
)

gene_truth <- data.frame(
  gene_id = gene_id,
  category = gene_category,
  nb_intercept = nb_parameters[, "intercept"],
  nb_group_2 = nb_parameters[, "group_2"],
  zi_intercept = zi_parameters[, "intercept"],
  zi_group_2 = zi_parameters[, "group_2"],
  dispersion = dispersion,
  stringsAsFactors = FALSE
)

saveRDS(simulated, file.path(output_dir, "simulated_dea.rds"))
utils::write.csv(
  simulated$object$count,
  file = file.path(output_dir, "count_matrix.csv")
)
utils::write.csv(
  cell_metadata,
  file = file.path(output_dir, "cell_metadata.csv"),
  row.names = FALSE
)
utils::write.csv(
  gene_truth,
  file = file.path(output_dir, "gene_truth.csv"),
  row.names = FALSE
)
utils::write.csv(simulated$mu, file = file.path(output_dir, "mu.csv"))
utils::write.csv(simulated$pi0, file = file.path(output_dir, "pi0.csv"))
utils::write.csv(
  simulated$structural_zero,
  file = file.path(output_dir, "structural_zero.csv")
)

cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat("Count dimensions:", paste(dim(simulated$object$count), collapse = " x "), "\n")
cat("Cell groups:", paste(names(table(group)), as.integer(table(group)), collapse = ", "), "\n")
cat(
  "Gene categories:",
  paste(names(table(gene_category)), as.integer(table(gene_category)), collapse = ", "),
  "\n"
)
