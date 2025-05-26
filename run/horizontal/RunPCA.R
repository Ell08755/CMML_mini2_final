# Seurat_horizontal_multi_batch_integration.R

# ----------- Load libraries -----------
library(Seurat)
library(Signac)
library(Matrix)
library(dplyr)
library(EnsDb.Hsapiens.v86)
library(harmony)

# ----------- Define input/output paths -----------
data_root <- "data/GSE156478/"
batches <- c("Control", "Stim")  # <-- modify this list
res_path <- "run_res/horizontal/GSE156478/"

# ----------- Initialize containers -----------
rna_list <- list()
atac_list <- list()

# ----------- Loop through batches -----------
for (batch in batches) {
  # --- Load RNA ---
  path_rna <- file.path(data_root, batch, "RNA")
  cell_names <- read.csv(file.path(path_rna, "barcodes.tsv"), header = FALSE)[[1]]
  gene_names <- read.table(file.path(path_rna, "features.tsv"), header = FALSE, sep = "\t")[[1]]
  counts_rna <- readMM(file.path(path_rna, "matrix.mtx"))
  rownames(counts_rna) <- gene_names
  colnames(counts_rna) <- cell_names

  rna <- CreateSeuratObject(counts = counts_rna, project = batch)
  rna$batch <- batch
  DefaultAssay(rna) <- "RNA"
  rna_list[[batch]] <- rna

  # --- Load ATAC ---
  path_atac <- file.path(data_root, batch, "ATAC")
  cell_names <- read.csv(file.path(path_atac, "barcodes.tsv"), header = FALSE)[[1]]
  peak_names <- read.table(file.path(path_atac, "features.tsv"), header = FALSE, sep = "\t")[[1]]
  counts_atac <- readMM(file.path(path_atac, "matrix.mtx"))
  rownames(counts_atac) <- peak_names
  colnames(counts_atac) <- cell_names

  chrom_assay <- CreateChromatinAssay(
    counts = counts_atac,
    sep = c("-", "-"),
    fragments = NULL
  )
  atac <- CreateSeuratObject(counts = counts_atac, assay = "ATAC", project = batch)
  atac[["ATAC"]] <- chrom_assay
  atac$batch <- batch
  DefaultAssay(atac) <- "ATAC"
  atac_list[[batch]] <- atac
}

# ----------- RNA Integration -----------
rna_list <- lapply(rna_list, function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, nfeatures = 10000)
})
features <- SelectIntegrationFeatures(rna_list, nfeatures = 4000)
rna_list <- lapply(rna_list, function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
})
anchors_rna <- FindIntegrationAnchors(object.list = rna_list, anchor.features = features, reduction = "rpca")
rna_merged <- IntegrateData(anchorset = anchors_rna)
rna_merged <- ScaleData(rna_merged)
rna_merged <- RunPCA(rna_merged)
rna_merged@assays

# ----------- ATAC Integration with Harmony -----------
atac_merged <- merge(atac_list[[1]], y = atac_list[-1])
atac_merged <- RunTFIDF(atac_merged)
atac_merged <- FindTopFeatures(atac_merged, min.cutoff = "q0")
atac_merged <- RunSVD(atac_merged)
# atac_merged <-RunHarmony(object=atac_merged, group.by.vars = "batch", reduction.use = "lsi", dims.use=2:50)
lsi=Embeddings(atac_merged, "lsi")[,2:50]
meta=atac_merged@meta.data
harmony_mat=HarmonyMatrix(
  data_mat = lsi,
  meta_data = meta,
  vars_use = "batch",
  do_pca=FALSE
)
atac_merged[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_mat,
  key="harmony_",
  assay="ATAC"
)


# ----------- Merge modalities and run joint PCA -----------
combined <- rna_merged
combined[["ATAC"]] <- atac_merged[["ATAC"]]
combined@reductions[["lsi"]] <- atac_merged@reductions[["harmony"]]

# Extract embeddings
rna_embed <- Embeddings(combined, reduction = "pca")[, 1:15]
atac_embed <- Embeddings(combined, reduction = "lsi")[, 1:15]
joint_embed <- cbind(rna_embed, atac_embed)
colnames(joint_embed)=paste0("Joint_",1:ncol(joint_embed))
# Create a new reduction object
combined[["combined"]] <- CreateDimReducObject(embeddings = joint_embed, key = "Joint_", assay = DefaultAssay(combined))

# Run PCA/UMAP/Clustering on combined
combined <- FindNeighbors(combined, reduction = "combined", dims = 1:30)
# combined <- FindClusters(combined, resolution = 0.5)
# combined <- RunUMAP(combined, reduction = "combined", dims = 1:30)

# ----------- Save output graph matrices -----------
writeMM(combined@graphs$integrated_snn, paste0(res_path,"PCA_connectivities.mtx"))
writeMM(combined@graphs$integrated_nn, paste0(res_path,"PCA_distance.mtx"))
write.csv(Embeddings(combined,"combined"), file=paste0(res_path,"PCA.csv"), quote=FALSE)
        