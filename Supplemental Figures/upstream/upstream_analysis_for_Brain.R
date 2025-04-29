### upstream analysis for Brain

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)
library(SoupX)
library(DropletUtils)
library(knitr)

### sample1
run_soupx <- function(toc, tod, rho = NULL){
  all <- toc
  all <- CreateSeuratObject(all)
  all <- NormalizeData(all, normalization.method = "LogNormalize", scale.factor = 10000)
  all <- FindVariableFeatures(all, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(all)
  all <- ScaleData(all, features = all.genes)
  all <- RunPCA(all, features = VariableFeatures(all), npcs = 30, verbose = F)
  all <- FindNeighbors(all, dims = 1:30)
  all <- FindClusters(all, resolution = 0.5)
  all <- RunUMAP(all, dims = 1:30)
  matx <- all@meta.data

  sc <- SoupChannel(tod, toc)
  sc <- setClusters(sc, setNames(matx$seurat_clusters, rownames(matx)))

  if (is.null(rho)){
   tryCatch({
     sc <- autoEstCont(sc)
   },
     error = function(e){
       print("autoEstCont Error !")
       sc <- setContaminationFraction(sc, 0.2)}
   ) } else {
   sc <- setContaminationFraction(sc, rho)
 }
  out <- adjustCounts(sc, roundToInt = TRUE)
  saveRDS(sc,"sample1.rds")
  DropletUtils:::write10xCounts("soupX_matrix", out, version = "3")
}


data_dir <- "RawMatrix"
list.files(data_dir)
tod <- Read10X(data.dir = data_dir, gene.column = 1)

data_dir <- "FilterMatrix"
list.files(data_dir)
toc <- Read10X(data.dir = data_dir, gene.column = 1)
tod <- tod[rownames(toc), ]
run_soupx(toc,tod, rho = 0.2)



brain_1 <- Read10X("soupX_matrix", gene.column = 1)
doublet <- read.table("sample1_Scrublet_result.txt", sep = "\t", header = T)
brain_1 <- brain_1[, !(colnames(brain_1) %in% doublet$CB)]

### Create Seurat Object and QC
brain_1 <- CreateSeuratObject(counts = brain_1, min.cells = 3, min.features = 200, project = "sample1")
brain_1[["percent.mt"]] <- PercentageFeatureSet(brain_1, pattern = "^mt-")

### Subset object according to nFeature_RNA, nCount_RNA, percent.mt
brain_1 <-subset(brain_1, subset = nFeature_RNA < 6000 & nFeature_RNA > 300 & percent.mt < 3)
brain_1@meta.data$orig.ident <- "sample1"