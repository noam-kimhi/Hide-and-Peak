
args <- commandArgs(trailingOnly = TRUE)

rds_path <- args[[1]]
output_dir <- args[[2]]

required_packages <- c("Seurat", "Signac")
missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste(
      "Missing R packages:",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before rerunning this script."
    )
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
})

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

message("Loading Seurat object...")

connection <- gzfile(rds_path, open = "rb")
object <- readRDS(connection)
close(connection)

message("Object loaded successfully.")

# Seurat's [[]] accessor returns the per-cell metadata.
metadata <- object[[]]

# Preserve identifiers explicitly as ordinary columns.
metadata$seurat_cell_id <- rownames(metadata)
metadata$active_ident <- as.character(Idents(object))

column_summary <- data.frame(
  column = colnames(metadata),
  r_class = vapply(
    metadata,
    function(column) paste(class(column), collapse = "/"),
    character(1)
  ),
  n_unique = vapply(
    metadata,
    function(column) length(unique(column)),
    integer(1)
  ),
  n_missing = vapply(
    metadata,
    function(column) sum(is.na(column)),
    integer(1)
  ),
  stringsAsFactors = FALSE
)

write.csv(
  column_summary,
  file = file.path(output_dir, "metadata_column_summary.csv"),
  row.names = FALSE
)

metadata_connection <- gzfile(
  file.path(output_dir, "GSE281367_cell_metadata.csv.gz"),
  open = "wt"
)

write.csv(
  metadata,
  file = metadata_connection,
  row.names = FALSE
)

close(metadata_connection)

candidate_columns <- grep(
  pattern = paste(
    "cell",
    "type",
    "annot",
    "cluster",
    "ident",
    "sample",
    "orig",
    "condition",
    "disease",
    "group",
    "sub",
    sep = "|"
  ),
  x = colnames(metadata),
  value = TRUE,
  ignore.case = TRUE
)

report <- capture.output({
  cat("=== SOFTWARE ===\n")
  cat("R version:", R.version.string, "\n")
  cat("Seurat version:", as.character(packageVersion("Seurat")), "\n")
  cat("Signac version:", as.character(packageVersion("Signac")), "\n\n")

  cat("=== OBJECT ===\n")
  cat("Class:", paste(class(object), collapse = ", "), "\n")
  cat("Object size:", format(object.size(object), units = "GB"), "\n")
  cat("Number of cells:", ncol(object), "\n")
  cat("Number of features in active assay:", nrow(object), "\n")
  cat("Default assay:", DefaultAssay(object), "\n")
  cat("Object slots:", paste(slotNames(object), collapse = ", "), "\n\n")

  cat("=== ASSAYS ===\n")
  assay_names <- Assays(object)
  cat("Assay names:", paste(assay_names, collapse = ", "), "\n")

  for (assay_name in assay_names) {
    assay_object <- object[[assay_name]]

    cat(
      "\nAssay:",
      assay_name,
      "-",
      nrow(assay_object),
      "features x",
      ncol(assay_object),
      "cells\n"
    )

    layers <- tryCatch(
      SeuratObject::Layers(assay_object),
      error = function(error) character(0)
    )

    if (length(layers) > 0) {
      cat("Layers:", paste(layers, collapse = ", "), "\n")
    }
  }

  cat("\n=== DIMENSIONAL REDUCTIONS ===\n")
  reduction_names <- Reductions(object)
  cat(
    "Reduction names:",
    paste(reduction_names, collapse = ", "),
    "\n"
  )

  for (reduction_name in reduction_names) {
    embedding <- Embeddings(object[[reduction_name]])

    cat(
      reduction_name,
      ":",
      nrow(embedding),
      "cells x",
      ncol(embedding),
      "dimensions\n"
    )
  }

  cat("\n=== METADATA ===\n")
  cat(
    "Metadata dimensions:",
    nrow(metadata),
    "rows x",
    ncol(metadata),
    "columns\n"
  )

  print(column_summary)

  cat("\n=== CANDIDATE ANNOTATION COLUMNS ===\n")
  print(candidate_columns)

  preview_columns <- unique(
    c("seurat_cell_id", "active_ident", candidate_columns)
  )

  cat("\n=== FIRST TEN CELLS ===\n")
  print(
    head(
      metadata[, preview_columns, drop = FALSE],
      10
    )
  )

  cat("\n=== VALUE COUNTS FOR CANDIDATE COLUMNS ===\n")

  for (column_name in candidate_columns) {
    cat("\n---", column_name, "---\n")

    value_counts <- sort(
      table(metadata[[column_name]], useNA = "ifany"),
      decreasing = TRUE
    )

    print(head(value_counts, 30))
  }

  cat("\n=== FIRST 30 SEURAT CELL IDENTIFIERS ===\n")
  print(head(metadata$seurat_cell_id, 30))
})

report_path <- file.path(
  output_dir,
  "seurat_object_report.txt"
)

writeLines(report, report_path)
cat(paste(report, collapse = "\n"))
cat("\n")
      reduction_name,
      ":",
      nrow(embedding),
      "cells x",
      ncol(embedding),
      "dimensions\n"
    )
  }

  cat("\n=== METADATA ===\n")
  cat(
    "Metadata dimensions:",
    nrow(metadata),
    "rows x",
    ncol(metadata),
    "columns\n"
  )

  print(column_summary)

  cat("\n=== CANDIDATE ANNOTATION COLUMNS ===\n")
  print(candidate_columns)

  preview_columns <- unique(
    c("seurat_cell_id", "active_ident", candidate_columns)
  )

  cat("\n=== FIRST TEN CELLS ===\n")
  print(
    head(
      metadata[, preview_columns, drop = FALSE],
      10
    )
  )

  cat("\n=== VALUE COUNTS FOR CANDIDATE COLUMNS ===\n")

  for (column_name in candidate_columns) {
    cat("\n---", column_name, "---\n")

    value_counts <- sort(
      table(metadata[[column_name]], useNA = "ifany"),
      decreasing = TRUE
    )

    print(head(value_counts, 30))
  }

  cat("\n=== FIRST 30 SEURAT CELL IDENTIFIERS ===\n")
  print(head(metadata$seurat_cell_id, 30))
})

report_path <- file.path(
  output_dir,
  "seurat_object_report.txt"
)

writeLines(report, report_path)
cat(paste(report, collapse = "\n"))
cat("\n")
