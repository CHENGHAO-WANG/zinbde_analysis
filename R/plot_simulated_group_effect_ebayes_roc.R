# Plot ROC curves for simulated group-effect Wald p-values before and after
# empirical Bayes dispersion shrinkage, using the known simulation truth.

required_packages <- c("ggplot2")
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

library(ggplot2)

wald_file <- file.path(
  "output",
  "simulated_group_effect_fit_settings_ebayes",
  "wald_p_values.csv"
)
truth_file <- file.path(
  "output",
  "simulated_group_effect_count_data",
  "gene_truth.csv"
)
output_dir <- file.path(
  "output",
  "simulated_group_effect_fit_settings_ebayes_roc"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stages <- c("before_ebayes", "after_ebayes")
components <- c("nb", "zi", "joint")
metric_specs <- data.frame(
  metric = c("p_value", "p_adj"),
  metric_label = c("nominal", "adjusted"),
  cutoff = c(0.05, 0.10),
  stringsAsFactors = FALSE
)

# Stop early if an input table is missing a column required by the analysis.
require_columns <- function(data, required_columns, table_name) {
  missing_columns <- setdiff(required_columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      sprintf(
        "%s is missing required column(s): %s.",
        table_name,
        paste(missing_columns, collapse = ", ")
      ),
      call. = FALSE
    )
  }
}

# Convert zinbde output gene labels such as gene1 to truth labels such as gene_001.
normalize_gene_id <- function(gene) {
  normalized <- as.character(gene)
  needs_padding <- grepl("^gene[0-9]+$", normalized)
  gene_number <- suppressWarnings(as.integer(sub("^gene", "", normalized[needs_padding])))
  normalized[needs_padding] <- sprintf("gene_%03d", gene_number)
  normalized
}

# Return the binary truth value corresponding to one Wald-test component.
component_truth <- function(component, truth) {
  if (identical(component, "nb")) {
    truth$nb_group_2 != 0
  } else if (identical(component, "zi")) {
    truth$zi_group_2 != 0
  } else if (identical(component, "joint")) {
    truth$nb_group_2 != 0 | truth$zi_group_2 != 0
  } else {
    stop(sprintf("Unsupported component: %s.", component), call. = FALSE)
  }
}

# Compute ROC points by ranking lower p-values as stronger evidence.
compute_roc_curve <- function(truth, p_value) {
  keep <- !is.na(truth) & is.finite(p_value) & p_value >= 0 & p_value <= 1
  truth <- as.logical(truth[keep])
  score <- 1 - p_value[keep]

  n_pos <- sum(truth)
  n_neg <- sum(!truth)
  if (n_pos == 0L || n_neg == 0L) {
    stop("ROC requires at least one positive and one negative gene.", call. = FALSE)
  }

  thresholds <- sort(unique(score), decreasing = TRUE)
  roc <- data.frame(
    threshold = c(Inf, thresholds),
    fpr = 0,
    tpr = 0
  )

  for (i in seq_along(thresholds)) {
    predicted <- score >= thresholds[i]
    roc$fpr[i + 1L] <- sum(predicted & !truth) / n_neg
    roc$tpr[i + 1L] <- sum(predicted & truth) / n_pos
  }

  roc
}

# Compute trapezoidal AUC from an ROC curve ordered by false-positive rate.
compute_auc <- function(roc) {
  roc <- roc[order(roc$fpr, roc$tpr), ]
  sum(diff(roc$fpr) * (head(roc$tpr, -1L) + tail(roc$tpr, -1L)) / 2)
}

# Summarize sensitivity and specificity at the requested p-value cutoff.
compute_cutoff_performance <- function(truth, p_value, cutoff) {
  keep <- !is.na(truth) & is.finite(p_value) & p_value >= 0 & p_value <= 1
  truth <- as.logical(truth[keep])
  predicted <- p_value[keep] <= cutoff

  tp <- sum(predicted & truth)
  fp <- sum(predicted & !truth)
  tn <- sum(!predicted & !truth)
  fn <- sum(!predicted & truth)

  data.frame(
    cutoff = cutoff,
    sensitivity = if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_,
    specificity = if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_,
    fpr = if ((fp + tn) > 0L) fp / (fp + tn) else NA_real_,
    tpr = if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_,
    true_positive = tp,
    false_positive = fp,
    true_negative = tn,
    false_negative = fn
  )
}

# Save one ROC plot for one stage, component, and p-value metric.
plot_roc <- function(roc_data, cutoff_data, stage, component, metric_label) {
  plot_data <- roc_data[
    roc_data$stage == stage &
      roc_data$component == component &
      roc_data$metric_label == metric_label,
    ,
    drop = FALSE
  ]
  point_data <- cutoff_data[
    cutoff_data$stage == stage &
      cutoff_data$component == component &
      cutoff_data$metric_label == metric_label,
    ,
    drop = FALSE
  ]

  roc_plot <- ggplot(
    plot_data,
    aes(x = fpr, y = tpr, color = setting, group = setting)
  ) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    geom_step(linewidth = 0.7) +
    geom_point(data = point_data, aes(x = fpr, y = tpr), size = 1.8) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    labs(
      x = "False positive rate",
      y = "True positive rate",
      color = "Fit setting",
      title = sprintf(
        "ROC: %s p-value, %s, %s",
        metric_label,
        gsub("_", " ", stage, fixed = TRUE),
        component
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold")
    )

  ggsave(
    filename = file.path(
      output_dir,
      sprintf("roc_%s_%s_%s.png", stage, component, metric_label)
    ),
    plot = roc_plot,
    width = 7,
    height = 5.5,
    dpi = 300
  )
}

cat("Wald p-value file:", normalizePath(wald_file, winslash = "/", mustWork = TRUE), "\n")
cat("Truth file:", normalizePath(truth_file, winslash = "/", mustWork = TRUE), "\n")

wald_results <- utils::read.csv(wald_file, stringsAsFactors = FALSE)
gene_truth <- utils::read.csv(truth_file, stringsAsFactors = FALSE)

require_columns(
  wald_results,
  c("setting", "stage", "gene", "component", "p_value", "p_adj"),
  "Wald p-value table"
)
require_columns(
  gene_truth,
  c("gene_id", "nb_group_2", "zi_group_2"),
  "Gene truth table"
)

wald_results$gene_id <- normalize_gene_id(wald_results$gene)
analysis_data <- merge(
  wald_results,
  gene_truth[, c("gene_id", "nb_group_2", "zi_group_2")],
  by = "gene_id",
  all.x = TRUE,
  sort = FALSE
)

if (anyNA(analysis_data$nb_group_2) || anyNA(analysis_data$zi_group_2)) {
  missing_gene <- unique(analysis_data$gene[is.na(analysis_data$nb_group_2)])
  stop(
    sprintf(
      "Could not match %d Wald gene(s) to the truth table. First missing gene: %s.",
      length(missing_gene),
      missing_gene[1]
    ),
    call. = FALSE
  )
}

analysis_data <- analysis_data[
  analysis_data$stage %in% stages & analysis_data$component %in% components,
  ,
  drop = FALSE
]

roc_rows <- list()
auc_rows <- list()
cutoff_rows <- list()
row_index <- 0L
auc_index <- 0L
cutoff_index <- 0L

settings <- sort(unique(analysis_data$setting))

for (stage in stages) {
  for (component in components) {
    for (metric_index in seq_len(nrow(metric_specs))) {
      metric <- metric_specs$metric[metric_index]
      metric_label <- metric_specs$metric_label[metric_index]
      cutoff <- metric_specs$cutoff[metric_index]

      for (setting in settings) {
        subset_data <- analysis_data[
          analysis_data$stage == stage &
            analysis_data$component == component &
            analysis_data$setting == setting,
          ,
          drop = FALSE
        ]

        truth <- component_truth(component, subset_data)
        roc <- compute_roc_curve(truth, subset_data[[metric]])
        auc <- compute_auc(roc)
        cutoff_performance <- compute_cutoff_performance(
          truth = truth,
          p_value = subset_data[[metric]],
          cutoff = cutoff
        )

        row_index <- row_index + 1L
        roc_rows[[row_index]] <- cbind(
          data.frame(
            setting = setting,
            stage = stage,
            component = component,
            metric = metric,
            metric_label = metric_label,
            stringsAsFactors = FALSE
          ),
          roc
        )

        auc_index <- auc_index + 1L
        auc_rows[[auc_index]] <- data.frame(
          setting = setting,
          stage = stage,
          component = component,
          metric = metric,
          metric_label = metric_label,
          auc = auc,
          stringsAsFactors = FALSE
        )

        cutoff_index <- cutoff_index + 1L
        cutoff_rows[[cutoff_index]] <- cbind(
          data.frame(
            setting = setting,
            stage = stage,
            component = component,
            metric = metric,
            metric_label = metric_label,
            stringsAsFactors = FALSE
          ),
          cutoff_performance
        )
      }

      plot_roc(
        roc_data = do.call(rbind, roc_rows),
        cutoff_data = do.call(rbind, cutoff_rows),
        stage = stage,
        component = component,
        metric_label = metric_label
      )
    }
  }
}

roc_summary <- do.call(rbind, roc_rows)
auc_summary <- do.call(rbind, auc_rows)
cutoff_performance <- do.call(rbind, cutoff_rows)

auc_summary$stage <- factor(auc_summary$stage, levels = stages)
auc_summary$component <- factor(auc_summary$component, levels = components)
auc_summary$metric_label <- factor(
  auc_summary$metric_label,
  levels = metric_specs$metric_label
)

auc_plot <- ggplot(
  auc_summary,
  aes(x = stage, y = auc, fill = setting)
) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_grid(component ~ metric_label) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    x = "Analysis stage",
    y = "AUC",
    fill = "Fit setting",
    title = "ROC AUC by empirical Bayes stage, component, and p-value type"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "auc_bar_plot.png"),
  plot = auc_plot,
  width = 10,
  height = 7,
  dpi = 300
)

utils::write.csv(
  roc_summary,
  file = file.path(output_dir, "roc_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  auc_summary,
  file = file.path(output_dir, "auc_summary.csv"),
  row.names = FALSE
)
utils::write.csv(
  cutoff_performance,
  file = file.path(output_dir, "cutoff_performance.csv"),
  row.names = FALSE
)

cat("Output directory:", normalizePath(output_dir, winslash = "/", mustWork = TRUE), "\n")
cat("ROC plots:", length(list.files(output_dir, pattern = "^roc_.*\\.png$")), "\n")
cat("AUC rows:", nrow(auc_summary), "\n")
