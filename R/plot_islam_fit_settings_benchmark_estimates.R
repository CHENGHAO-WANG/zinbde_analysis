# Plot coefficient, SE, and dispersion agreement across Islam benchmark
# settings. The script excludes genes with discordant zero-inflation decisions,
# keeps consistently ZINB genes, and writes centered-difference and CV plots.

benchmark_dir <- file.path("output", "islam_fit_settings_benchmark")
plot_dir <- file.path("output", "islam_fit_settings_benchmark_estimate_plots")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

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

# Generate enough distinct colors for setting and estimate-series lines.
line_colors <- function(n) {
  grDevices::hcl.colors(n, palette = "Dark 3")
}

# Draw a legend in its own plotting panel so long labels are not clipped by the
# main plot margins or the PNG device boundary.
draw_legend_panel <- function(labels, colors, columns) {
  par(mar = c(0, 5, 0, 2), xpd = NA)
  plot.new()
  legend(
    "center",
    legend = labels,
    col = colors,
    lty = 1,
    lwd = 1.2,
    ncol = columns,
    cex = 0.58,
    bty = "n"
  )
}

# Plot centered differences for one estimate series. Genes are ordered by that
# series' CV so each plot uses its own stability ranking.
plot_centered_difference <- function(series_name) {
  series_data <- estimate_data[estimate_data$series == series_name, ]
  gene_order <- series_cv_order(series_name)
  series_data$gene_index <- match(series_data$gene_id, gene_order)
  series_data <- series_data[order(series_data$setting, series_data$gene_index), ]

  output_file <- file.path(
    plot_dir,
    sprintf("centered_difference_%s.png", safe_filename(series_name))
  )

  grDevices::png(output_file, width = 2800, height = 1700, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1.3))
  par(mar = c(8, 5, 4, 2), xpd = FALSE)
  y_range <- range(series_data$difference, na.rm = TRUE)

  # Initialize the plot without data so all settings share the same axes.
  plot(
    seq_along(gene_order),
    rep(NA_real_, length(gene_order)),
    type = "n",
    xaxt = "n",
    xlab = "Genes ordered by CV for this estimate series",
    ylab = "Estimate minus gene-wise setting mean",
    ylim = y_range,
    main = sprintf("Centered differences: %s", series_name)
  )
  axis(
    1,
    at = seq_along(gene_order),
    labels = gene_order,
    las = 2,
    cex.axis = 0.65
  )
  abline(h = 0, col = "gray70", lty = 2)

  # Draw one line per benchmark setting using the series-specific gene order.
  plot_settings <- sort(unique(series_data$setting))
  colors <- line_colors(length(plot_settings))
  for (i in seq_along(plot_settings)) {
    setting_data <- series_data[series_data$setting == plot_settings[i], ]
    setting_data <- setting_data[order(setting_data$gene_index), ]
    lines(
      setting_data$gene_index,
      setting_data$difference,
      col = colors[i],
      lwd = 1.2
    )
  }
  draw_legend_panel(plot_settings, colors, columns = 2L)

  invisible(output_file)
}

# Plot CV for all estimate series. Genes are ordered by their average CV across
# available coefficient, SE, and dispersion series.
plot_cv <- function() {
  cv_data <- summary_values[summary_values$gene_id %in% average_cv_order, ]
  cv_data$gene_index <- match(cv_data$gene_id, average_cv_order)
  cv_data <- cv_data[order(cv_data$series, cv_data$gene_index), ]
  series_names <- sort(unique(cv_data$series))

  output_file <- file.path(plot_dir, "coefficient_of_variation_by_gene.png")

  grDevices::png(output_file, width = 2800, height = 1700, res = 180)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  })

  layout(matrix(c(1, 2), nrow = 2), heights = c(4, 1.1))
  par(mar = c(8, 5, 4, 2), xpd = FALSE)
  y_range <- range(cv_data$cv, na.rm = TRUE)

  # Initialize the plot without data so all estimate series share the same axes.
  plot(
    seq_along(average_cv_order),
    rep(NA_real_, length(average_cv_order)),
    type = "n",
    xaxt = "n",
    xlab = "Genes ordered by average CV",
    ylab = "Coefficient of variation",
    ylim = y_range,
    main = "Coefficient of variation across settings"
  )
  axis(
    1,
    at = seq_along(average_cv_order),
    labels = average_cv_order,
    las = 2,
    cex.axis = 0.65
  )

  # Draw one line per estimate series using the average-CV gene order.
  colors <- line_colors(length(series_names))
  for (i in seq_along(series_names)) {
    series_data <- cv_data[cv_data$series == series_names[i], ]
    series_data <- series_data[order(series_data$gene_index), ]
    lines(
      series_data$gene_index,
      series_data$cv,
      col = colors[i],
      lwd = 1.2
    )
  }
  draw_legend_panel(series_names, colors, columns = 3L)

  invisible(output_file)
}

plot_files <- c(
  vapply(sort(unique(estimate_data$series)), plot_centered_difference, character(1)),
  plot_cv()
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
