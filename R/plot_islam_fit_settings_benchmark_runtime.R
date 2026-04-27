# Plot elapsed benchmark runtime by setting. The script reads the runtime
# summary written by islam_fit_settings_benchmark.R, orders settings from
# fastest to slowest, and writes a PNG bar plot beside the runtime CSV.

required_packages <- c("ggplot2")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    sprintf(
      "Install required package(s) before plotting runtime: %s.",
      paste(missing_packages, collapse = ", ")
    ),
    call. = FALSE
  )
}

benchmark_dir <- file.path("output", "islam_fit_settings_benchmark")
runtime_file <- file.path(benchmark_dir, "runtime.csv")
plot_file <- file.path(benchmark_dir, "runtime_elapsed_seconds.png")

if (!file.exists(runtime_file)) {
  stop(sprintf("Missing runtime file: %s.", runtime_file), call. = FALSE)
}

runtime <- read.csv(runtime_file, stringsAsFactors = FALSE)
required_columns <- c("setting", "elapsed_seconds")
missing_columns <- setdiff(required_columns, names(runtime))
if (length(missing_columns) > 0L) {
  stop(
    sprintf(
      "runtime.csv is missing required column(s): %s.",
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

runtime$elapsed_seconds <- as.numeric(runtime$elapsed_seconds)
runtime <- runtime[order(runtime$elapsed_seconds), ]
runtime$setting <- factor(runtime$setting, levels = runtime$setting)

runtime_plot <- ggplot2::ggplot(
  runtime,
  ggplot2::aes(x = setting, y = elapsed_seconds)
) +
  ggplot2::geom_col(fill = "#3B82F6", width = 0.75) +
  ggplot2::coord_flip() +
  ggplot2::labs(
    x = "Setting",
    y = "Elapsed seconds",
    title = "Runtime by benchmark setting"
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid.major.y = ggplot2::element_blank()
  )

ggplot2::ggsave(
  filename = plot_file,
  plot = runtime_plot,
  width = 10,
  height = max(4, 0.25 * nrow(runtime)),
  dpi = 300
)
