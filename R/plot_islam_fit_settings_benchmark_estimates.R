# Plot coefficient, SE, and dispersion agreement across Islam benchmark
# settings. The script excludes genes with discordant zero-inflation decisions,
# keeps consistently ZINB genes, and writes centered-difference and CV heatmaps.

benchmark_dir <- file.path("output", "islam_fit_settings_benchmark")
plot_dir <- file.path("output", "islam_fit_settings_benchmark_estimate_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(plot_dir, pattern = "[.]png$", full.names = TRUE))

discordant_file <- file.path(benchmark_dir, "discordant_zi_decision_model_types.csv")
decisions_file <- file.path(benchmark_dir, "zi_decisions.csv")
coefficients_file <- file.path(benchmark_dir, "coefficients.csv")
dispersion_file <- file.path(benchmark_dir, "dispersion.csv")

required_files <- c(
  discordant_file,
  decisions_file,
  coefficients_file,
  dispersion_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop(
    sprintf("Missing required benchmark file(s): %s.", paste(missing_files, collapse = ", ")),
    call. = FALSE
  )
}

discordant_model_types <- read.csv(discordant_file, check.names = FALSE)
zi_decisions <- read.csv(decisions_file, stringsAsFactors = FALSE)
coefficients <- read.csv(coefficients_file, stringsAsFactors = FALSE)
dispersion <- read.csv(dispersion_file, stringsAsFactors = FALSE)

discordant_genes <- setdiff(names(discordant_model_types), "setting")
filtered_decisions <- zi_decisions[!zi_decisions$gene_id %in% discordant_genes, ]
settings <- sort(unique(filtered_decisions$setting))

decision_counts <- xtabs(~ gene_id + model_type, data = filtered_decisions)
if (!"ZINB" %in% colnames(decision_counts)) {
  stop("No ZINB model decisions were found after filtering discordant genes.", call. = FALSE)
}

zinb_genes <- rownames(decision_counts)[
  decision_counts[, "ZINB"] == length(settings) &
    rowSums(decision_counts) == length(settings)
]
if (length(zinb_genes) == 0L) {
  stop("No consistently ZINB genes remained after filtering.", call. = FALSE)
}

# Build stable labels that identify a plotted estimate series.
make_series_id <- function(prefix, component = NULL, term = NULL) {
  values <- c(prefix, component, term)
  values <- values[!is.null(values) & !is.na(values)]
  paste(values, collapse = "__")
}

# Convert estimate series labels into portable plot filenames.
safe_filename <- function(value) {
  value <- gsub("[^A-Za-z0-9]+", "_", value)
  gsub("(^_+|_+$)", "", value)
}

# Convert coefficient or SE rows into the common long estimate format.
coefficient_metric_frame <- function(data, value_col, metric_name) {
  data.frame(
    setting = data$setting,
    gene_id = data$gene_id,
    series = mapply(
      make_series_id,
      prefix = metric_name,
      component = data$component,
      term = data$term,
      USE.NAMES = FALSE
    ),
    value = data[[value_col]],
    stringsAsFactors = FALSE
  )
}

filtered_coefficients <- coefficients[coefficients$gene_id %in% zinb_genes, ]
filtered_dispersion <- dispersion[dispersion$gene_id %in% zinb_genes, ]

estimate_data <- rbind(
  coefficient_metric_frame(filtered_coefficients, "estimate", "coefficient"),
  coefficient_metric_frame(filtered_coefficients, "std_error", "se"),
  data.frame(
    setting = filtered_dispersion$setting,
    gene_id = filtered_dispersion$gene_id,
    series = "dispersion",
    value = filtered_dispersion$dispersion,
    stringsAsFactors = FALSE
  )
)

estimate_data <- estimate_data[complete.cases(estimate_data[c("setting", "gene_id", "series")]), ]
estimate_data$value <- as.numeric(estimate_data$value)

# Return the across-setting mean and CV for a gene-estimate series.
compute_summary <- function(values) {
  estimate_mean <- mean(values, na.rm = TRUE)
  estimate_sd <- stats::sd(values, na.rm = TRUE)
  non_missing <- sum(!is.na(values))
  cv <- if (non_missing > 1L && is.finite(estimate_mean) && estimate_mean != 0) {
    estimate_sd / abs(estimate_mean)
  } else {
    NA_real_
  }

  c(mean = estimate_mean, cv = cv)
}

summary_values <- aggregate(
  value ~ gene_id + series,
  data = estimate_data,
  FUN = compute_summary
)
summary_values$mean <- summary_values$value[, "mean"]
summary_values$cv <- summary_values$value[, "cv"]
summary_values$value <- NULL

estimate_data <- merge(
  estimate_data,
  summary_values,
  by = c("gene_id", "series"),
  all.x = TRUE,
  sort = FALSE
)
estimate_data$difference <- estimate_data$value - estimate_data$mean

average_cv <- aggregate(cv ~ gene_id, data = summary_values, FUN = mean, na.rm = TRUE)
average_cv <- average_cv[order(average_cv$cv, average_cv$gene_id), ]
average_cv_order <- average_cv$gene_id

# Order genes by the CV for one centered-difference estimate series.
series_cv_order <- function(series_name) {
  series_summary <- summary_values[summary_values$series == series_name, ]
  series_summary <- series_summary[order(series_summary$cv, series_summary$gene_id, na.last = TRUE), ]
  series_summary$gene_id
}

# Build a numeric matrix from long-form data using explicit row and column
# orders. Duplicate row-column pairs are not expected in the benchmark summaries.
long_to_matrix <- function(data, row_col, col_col, value_col, row_order, col_order) {
  heatmap_matrix <- matrix(
    NA_real_,
    nrow = length(row_order),
    ncol = length(col_order),
    dimnames = list(row_order, col_order)
  )
  row_index <- match(data[[row_col]], row_order)
  col_index <- match(data[[col_col]], col_order)
  keep <- !is.na(row_index) & !is.na(col_index)
  heatmap_matrix[cbind(row_index[keep], col_index[keep])] <- data[[value_col]][keep]

  heatmap_matrix
}

# Cluster heatmap rows while preserving the user-defined gene order on columns.
cluster_row_order <- function(heatmap_matrix) {
  if (nrow(heatmap_matrix) < 2L) {
    return(seq_len(nrow(heatmap_matrix)))
  }

  clustering_matrix <- heatmap_matrix
  for (row in seq_len(nrow(clustering_matrix))) {
    missing <- is.na(clustering_matrix[row, ])
    row_mean <- mean(clustering_matrix[row, ], na.rm = TRUE)
    if (!is.finite(row_mean)) {
      row_mean <- 0
    }
    clustering_matrix[row, missing] <- row_mean
  }

  stats::hclust(stats::dist(clustering_matrix))$order
}

# Make a symmetric blue-white-red scale for centered differences.
difference_breaks <- function(values, color_count) {
  max_abs <- max(abs(values), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) {
    max_abs <- 1
  }
  seq(-max_abs, max_abs, length.out = color_count + 1L)
}

# Make a white-yellow-red scale for non-negative CV values.
cv_breaks <- function(values, color_count) {
  max_value <- max(values, na.rm = TRUE)
  if (!is.finite(max_value) || max_value == 0) {
    max_value <- 1
  }
  seq(0, max_value, length.out = color_count + 1L)
}

# Draw the color scale below each heatmap so the data panel can keep full labels.
draw_color_key <- function(colors, breaks, label) {
  key_values <- seq(min(breaks), max(breaks), length.out = length(colors))

  par(mar = c(5, 28, 0.5, 4))
  image(
    x = key_values,
    y = 1,
    z = matrix(key_values, nrow = length(key_values), ncol = 1L),
    col = colors,
    breaks = breaks,
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(1, cex.axis = 0.75)
  mtext(label, side = 1, line = 2.7, cex = 0.8)
}

# Draw a heatmap with clustered rows and fixed gene columns.
plot_heatmap <- function(
    heatmap_matrix,
    colors,
    breaks,
    output_file,
    main,
    xlab,
    key_label,
    y_axis_cex = 0.5) {
  row_order <- cluster_row_order(heatmap_matrix)
  heatmap_matrix <- heatmap_matrix[row_order, , drop = FALSE]

  grDevices::png(output_file, width = 3600, height = 2200, res = 220)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  layout(matrix(c(1, 2), nrow = 2), heights = c(6, 0.8))
  par(mar = c(8, 28, 4, 4))

  # Reverse rows for image() so the first clustered row appears at the top.
  display_matrix <- heatmap_matrix[rev(seq_len(nrow(heatmap_matrix))), , drop = FALSE]
  image(
    x = seq_len(ncol(display_matrix)),
    y = seq_len(nrow(display_matrix)),
    z = t(display_matrix),
    col = colors,
    breaks = breaks,
    axes = FALSE,
    xlab = xlab,
    ylab = "",
    main = main
  )
  axis(
    1,
    at = seq_len(ncol(display_matrix)),
    labels = colnames(display_matrix),
    las = 2,
    cex.axis = 0.62
  )
  axis(
    2,
    at = seq_len(nrow(display_matrix)),
    labels = rownames(display_matrix),
    las = 1,
    cex.axis = y_axis_cex
  )
  box()

  draw_color_key(colors, breaks, key_label)

  invisible(output_file)
}

# Plot centered differences for one estimate series. Genes are ordered by that
# series' CV, while settings on the other axis are clustered.
plot_centered_difference_heatmap <- function(series_name) {
  series_data <- estimate_data[estimate_data$series == series_name, ]
  gene_order <- series_cv_order(series_name)
  plot_settings <- sort(unique(series_data$setting))
  heatmap_matrix <- long_to_matrix(
    series_data,
    row_col = "setting",
    col_col = "gene_id",
    value_col = "difference",
    row_order = plot_settings,
    col_order = gene_order
  )

  output_file <- file.path(
    plot_dir,
    sprintf("heatmap_centered_difference_%s.png", safe_filename(series_name))
  )
  colors <- grDevices::colorRampPalette(c("#2166ac", "#f7f7f7", "#b2182b"))(101)

  plot_heatmap(
    heatmap_matrix = heatmap_matrix,
    colors = colors,
    breaks = difference_breaks(series_data$difference, length(colors)),
    output_file = output_file,
    main = sprintf("Centered differences: %s", series_name),
    xlab = "Genes ordered by CV for this estimate series",
    key_label = "Estimate minus gene-wise setting mean",
    y_axis_cex = 0.42
  )
}

# Plot CV for all estimate series. Genes are ordered by average CV, while
# estimate series on the other axis are clustered.
plot_cv_heatmap <- function() {
  cv_data <- summary_values[summary_values$gene_id %in% average_cv_order, ]
  series_names <- sort(unique(cv_data$series))
  heatmap_matrix <- long_to_matrix(
    cv_data,
    row_col = "series",
    col_col = "gene_id",
    value_col = "cv",
    row_order = series_names,
    col_order = average_cv_order
  )
  output_file <- file.path(plot_dir, "heatmap_coefficient_of_variation_by_gene.png")
  colors <- grDevices::colorRampPalette(c("#fffff7", "#fed976", "#bd0026"))(101)

  plot_heatmap(
    heatmap_matrix = heatmap_matrix,
    colors = colors,
    breaks = cv_breaks(cv_data$cv, length(colors)),
    output_file = output_file,
    main = "Coefficient of variation across settings",
    xlab = "Genes ordered by average CV",
    key_label = "Coefficient of variation",
    y_axis_cex = 0.58
  )
}

plot_files <- c(
  vapply(sort(unique(estimate_data$series)), plot_centered_difference_heatmap, character(1)),
  plot_cv_heatmap()
)

write.csv(
  summary_values[order(summary_values$series, summary_values$gene_id), ],
  file.path(plot_dir, "estimate_cv_summary.csv"),
  row.names = FALSE
)
write.csv(
  data.frame(
    gene_id = average_cv_order,
    average_cv = average_cv$cv,
    stringsAsFactors = FALSE
  ),
  file.path(plot_dir, "average_cv_gene_order.csv"),
  row.names = FALSE
)

message(sprintf("Excluded %d discordant gene(s).", length(discordant_genes)))
message(sprintf("Plotted %d consistently ZINB gene(s).", length(zinb_genes)))
message(sprintf("Wrote %d plot file(s) to %s.", length(plot_files), plot_dir))
