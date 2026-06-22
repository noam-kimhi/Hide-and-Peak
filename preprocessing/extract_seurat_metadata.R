library(Seurat)

rds_path <- "data/snatac_gse281367/GSE281367_seurat_clustered.rds.gz"
out_path  <- "data/derived/seurat_cell_metadata.csv"

cat("Loading Seurat object from", rds_path, "...\n")
obj <- readRDS(gzcon(file(rds_path, "rb")))

cat("Object class:", class(obj), "\n")
cat("Number of cells:", ncol(obj), "\n")
cat("Number of features:", nrow(obj), "\n\n")

meta <- obj@meta.data
cat("Metadata columns:\n")
print(colnames(meta))
cat("\nFirst 5 rows:\n")
print(head(meta, 5))

write.csv(meta, out_path, quote = FALSE)
cat("\nMetadata saved to", out_path, "\n")
