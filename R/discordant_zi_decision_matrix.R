# Summarize zero-inflation decision disagreements across benchmark settings.
# The output keeps settings on rows, discordant genes on columns, and uses the
# selected model type (NB or ZINB) as each cell value.

input_file <- file.path("output", "islam_fit_settings_benchmark", "zi_decisions.csv")
output_file <- file.path(
  "output",
  "islam_fit_settings_benchmark",
  "discordant_zi_decision_model_types.csv"
)

zi_decisions <- utils::read.csv(input_file, stringsAsFactors = FALSE)

required_columns <- c("setting", "gene_id", "model_type", "needs_zero_inflation")
missing_columns <- setdiff(required_columns, names(zi_decisions))
if (length(missing_columns) > 0L) {
  stop(
    sprintf(
      "Missing required column(s) in %s: %s.",
      input_file,
      paste(missing_columns, collapse = ", ")
    ),
    call. = FALSE
  )
}

decision_by_gene <- split(zi_decisions$needs_zero_inflation, zi_decisions$gene_id)
discordant_genes <- names(decision_by_gene)[
  vapply(decision_by_gene, function(x) length(unique(x)) > 1L, logical(1))
]

discordant_rows <- zi_decisions[zi_decisions$gene_id %in% discordant_genes, ]

setting_levels <- unique(zi_decisions$setting)
gene_levels <- unique(zi_decisions$gene_id[zi_decisions$gene_id %in% discordant_genes])

model_type_matrix <- tapply(
  discordant_rows$model_type,
  list(
    setting = factor(discordant_rows$setting, levels = setting_levels),
    gene_id = factor(discordant_rows$gene_id, levels = gene_levels)
  ),
  function(x) paste(unique(x), collapse = ";")
)

discordant_model_types <- data.frame(
  setting = rownames(model_type_matrix),
  model_type_matrix,
  row.names = NULL,
  check.names = FALSE
)

utils::write.csv(discordant_model_types, output_file, row.names = FALSE)

cat("Discordant gene count:", length(discordant_genes), "\n")
cat("Output file:", normalizePath(output_file, winslash = "/", mustWork = TRUE), "\n")
