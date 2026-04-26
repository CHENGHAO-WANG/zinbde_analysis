#######################################
# Load expression data

# Read expression .txt file from the data directory.
data_dir <- file.path(getwd(), "data/raw/GSE29087_L139_expression_tab.txt")

# The top lines are metadata (including Barcode and Sample rows).
metadata_lines <- readLines(data_dir, n = 6, warn = FALSE)
sample_tokens <- strsplit(metadata_lines[5], "\t", fixed = TRUE)[[1]]
sample_start <- match("Sample:", sample_tokens)
sample_ids <- sample_tokens[(sample_start + 1):length(sample_tokens)]

# Main data starts after the 6 metadata/header lines.
expression_data <- read.delim(
  data_dir,
  sep = "\t",
  header = FALSE,
  skip = 6,
  fill = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

base_cols <- strsplit(metadata_lines[6], "\t", fixed = TRUE)[[1]]
  
colnames(expression_data) <- c(base_cols, sample_ids)

rownames(expression_data) <- expression_data$Feature

expression_data[, base_cols] <- list(NULL)

# view the data
expression_data[1:10, 1:10]

#######################################
# Create metadata

meta_data <- data.frame(cell_id = colnames(expression_data))

# According to https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE29087
# Samples section
# Embryonic stem cell: samples 1-48 (A01 to H06), GEO GSM729107-GSM729154 (L139_ESC_1-L139_ESC_48)
# Embryonic fibroblast: samples 49-92 (A07 to D12), GEO GSM729155-GSM729198 (L139_MEF_49-L139_MEF_92)
# Negative control: samples 93-96 (E12, F12, G12, H12), GEO GSM729199-GSM729202 (L139_EMPTY_93-L139_EMPTY_96)

meta_data$label <- ""
meta_data$label[1:48] <- "Embryonic stem cell"
meta_data$label[49:92] <- "Embryonic fibroblast"
meta_data$label[93:96] <- "Negative control"

rownames(meta_data) <- meta_data$cell_id


library(zinbde)

islam <- create_dea_object(expression_data, meta_data)

dim(islam)

# remove negative control from cells
islam <- subset(islam, subset = label != "Negative control")
# remove RNA spikes from genes
islam <- islam[!grepl("^RNA_SPIKE", rownames(islam$count)),]

dim(islam)

islam[["cell_id"]] <- NULL

head(islam$metadata)

head(rownames(islam))

# save
saveRDS(islam, file = "data/processed/islam_dea_obj.rds")

