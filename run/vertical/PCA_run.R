# better run in powershell
# argument 1 is the path to the input file family, no / end
# argument 2 save path of the results, end will automatically add MIRA.csv

library(Seurat)
require(Matrix)
library(dplyr)
require(stringr)
library(Signac)
library(EnsDb.Hsapiens.v86)

args <- commandArgs(trailingOnly = TRUE)
# args=c("data/brain_SNARE", "run_res/vertical/brain_SNARE/")

# check existence
if (!file.exists(args[1])) {
  stop("Input data does not exist.")
}

cat("Input file:", args[1], "\n")

### RNA ###
path = paste0(args[1],'/RNA/')
Cell_name <- read.csv(paste0(path,'barcodes.tsv'),header = F)
Gene_name <- read.table(paste0(path,"features.tsv"),header = F, sep = "\t")
M <- readMM(paste0(path,"matrix.mtx"))
rownames(M) = Gene_name$V1
colnames(M) = Cell_name$V1

pbmc <- CreateSeuratObject(counts = M, project = "DOGMA", min.cells = 1, min.features = 1)
rm(M)
DefaultAssay(pbmc) <- "RNA"
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, nfeatures = 4000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, npcs = 50, reduction.name = "pca")

### ATAC ###
path = paste0(args[1],'/ATAC/')
Cell_name <- read.csv(paste0(path,'barcodes.tsv'),header = F)
Gene_name <- read.table(paste0(path,"features.tsv"),header = F, sep = "\t")
atac_counts <- readMM(paste0(path,"matrix.mtx"))
rownames(atac_counts) = Gene_name$V1
colnames(atac_counts) = Cell_name$V1

chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c("-", "-"),
  fragments = NULL
)
pbmc[["ATAC"]] <- chrom_assay
DefaultAssay(pbmc) <- "ATAC"
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc, n = 50, reduction.name = "lsi")

### concat pca and lsi ###
rna_embed <- Embeddings(pbmc, reduction = "pca")[,1:15]
atac_embed <- Embeddings(pbmc, reduction = "lsi")[,1:15]
joint_embed <- cbind(rna_embed, atac_embed)
joint_embed <- joint_embed[colnames(pbmc),,drop=FALSE]
joint_embed <- joint_embed[order(rownames(joint_embed)),,drop=FALSE]
pbmc=pbmc[,rownames(joint_embed)]
colnames(joint_embed)=paste0("Joint_",1:ncol(joint_embed))
pbmc[["combined"]] <- CreateDimReducObject(embeddings = joint_embed, key = "Joint_", assay = DefaultAssay(pbmc))

### combined ###
pbmc <- FindNeighbors(pbmc, reduction = "combined", dims = 1:ncol(joint_embed))

# pbmc <- FindClusters(pbmc, resolution = 0.5)
# pbmc <- RunUMAP(pbmc, reduction = "combined", dims = 1:ncol(joint_embed), reduction.name = "pca_lsi.umap", reduction.key = "pcaLSI_")

### save ###
dir <- strsplit(args[1],'/')[[1]]
writeMM(pbmc@graphs$ATAC_snn, paste0(args[2],"PCA_connectivities.mtx"))
writeMM(pbmc@graphs$ATAC_nn, paste0(args[2],"PCA_distance.mtx"))
