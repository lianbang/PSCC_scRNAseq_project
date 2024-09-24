
library(SeuratObject) # v4.1.3
library(Seurat) # v4.3.0
library(MySeuratWrappers)
library(ggplot2)
library(paletteer)

### Rscirpt for Figure S2A, S2C

### read cellranger output files
# All genes: 36601
P1N.data <- Read10X(data.dir = "cellranger-6.0.1/P1N-10X3/filtered_feature_bc_matrix/") # 8836 raw cells
P1T.data <- Read10X(data.dir = "cellranger-6.0.1/P1T-10X3/filtered_feature_bc_matrix/") # 10009 raw cells
P2N.data <- Read10X(data.dir = "cellranger-6.0.1/P2N-10X3/filtered_feature_bc_matrix/") # 4260 raw cells
P2T.data <- Read10X(data.dir = "cellranger-6.0.1/P2T-10X3/filtered_feature_bc_matrix/") # 5573 raw cells
P3N.data <- Read10X(data.dir = "cellranger-6.0.1/P3N-10X3/filtered_feature_bc_matrix/") # 3002 raw cells
P3T.data <- Read10X(data.dir = "cellranger-6.0.1/P3T-10X3/filtered_feature_bc_matrix/") # 5995 raw cells
P4T.data <- Read10X(data.dir = "cellranger-6.0.1/P4T-10X3/filtered_feature_bc_matrix/") # 4571 raw cells
P5N.data <- Read10X(data.dir = "cellranger-6.0.1/P5N-10X3/filtered_feature_bc_matrix/") # 4191 raw cells
P5T.data <- Read10X(data.dir = "cellranger-6.0.1/P5T-10X3/filtered_feature_bc_matrix/") # 5687 raw cells
P6N.data <- Read10X(data.dir = "cellranger-6.0.1/P6N-10X3/filtered_feature_bc_matrix/") # 11269 raw cells
P6T.data <- Read10X(data.dir = "cellranger-6.0.1/P6T-10X3/filtered_feature_bc_matrix/") # 8305 raw cells
PT16P.data <- Read10X(data.dir = "cellranger-6.0.1/PT16P-10X5/filtered_feature_bc_matrix/") # 19026 raw cells
PT19P.data <- Read10X(data.dir = "cellranger-6.0.1/PT19P-10X5/filtered_feature_bc_matrix/") # 21308 raw cells
PT2N.data <- Read10X(data.dir = "cellranger-6.0.1/PT2N-10X5/filtered_feature_bc_matrix/") # 18147 raw cells
PT6N.data <- Read10X(data.dir = "cellranger-6.0.1/PT6N-10X5/filtered_feature_bc_matrix/") # 19517 raw cells

### create seurat object
P1N <- CreateSeuratObject(counts = P1N.data, min.cells = 10, min.features = 200, project = "PN1") # 21247 features, 8778 cells
P1T <- CreateSeuratObject(counts = P1T.data, min.cells = 10, min.features = 200, project = "PT1") # 19484 features, 9969 cells
P2N <- CreateSeuratObject(counts = P2N.data, min.cells = 10, min.features = 200, project = "PN2") # 19395 features, 4237 cells
P2T <- CreateSeuratObject(counts = P2T.data, min.cells = 10, min.features = 200, project = "PT2") # 19612 features, 5443 cells
P3N <- CreateSeuratObject(counts = P3N.data, min.cells = 10, min.features = 200, project = "PN3") # 18647 features, 2996 cells
P3T <- CreateSeuratObject(counts = P3T.data, min.cells = 10, min.features = 200, project = "PT3") # 19460 features, 5978 cells
P4T <- CreateSeuratObject(counts = P4T.data, min.cells = 10, min.features = 200, project = "PT4") # 19479 features, 4292 cells
P5N <- CreateSeuratObject(counts = P5N.data, min.cells = 10, min.features = 200, project = "PN5") # 18518 features, 3974 cells
P5T <- CreateSeuratObject(counts = P5T.data, min.cells = 10, min.features = 200, project = "PT5") # 19577 features, 5575 cells
P6N <- CreateSeuratObject(counts = P6N.data, min.cells = 10, min.features = 200, project = "PN6") # 21322 features, 11239 cells
P6T <- CreateSeuratObject(counts = P6T.data, min.cells = 10, min.features = 200, project = "PT6") # 20533 features, 8291 cells
PT6N <- CreateSeuratObject(counts = PT6N.data, min.cells = 10, min.features = 200, project = "PT6N") # 16492 features, 19513 cells
PT2N <- CreateSeuratObject(counts = PT2N.data, min.cells = 10, min.features = 200, project = "PT2N")  # 16508 features, 18121 cells
PT19P <- CreateSeuratObject(counts = PT19P.data, min.cells = 10, min.features = 200, project = "PT19P")  # 17009 features, 21296 cells
PT16P <- CreateSeuratObject(counts = PT16P.data, min.cells = 10, min.features = 200, project = "PT16P")  # 17500 features, 18994 cells

PSCC.raw <- merge(x = P1N, y = c(P1T,P2N,P2T,P3N,P3T,P4T,P5N,P5T,P6N,P6T,PT2N,PT6N,PT16P,PT19P),add.cell.ids = c("P1N","P1T","P2N","P2T","P3N","P3T","P4T","P5N","P5T","P6N","P6T","PT2N","PT6N","PT16P","PT19P")) # 25312 features, 148696 cells
PSCC.raw[["percent.mt"]] <- PercentageFeatureSet(PSCC.raw, pattern = "MT-")
PSCC  <- subset(PSCC.raw, subset = nFeature_RNA > 500 & nFeature_RNA < 5000 & percent.mt < 10)  # 25312 features, 119351 cells

### add info into meta.data
PSCC@meta.data$Patient <- "Patient"
PSCC@meta.data$Tissue <- "Tissue"
PSCC@meta.data$HPV <- "HPV"
PSCC@meta.data$Groups <- "Groups"
PSCC@meta.data$Library <- "Library"

PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PN1" | orig.ident=="PT1")),]$Patient <- "P1"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PN2" | orig.ident=="PT2")),]$Patient <- "P2"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PN3" | orig.ident=="PT3")),]$Patient <- "P3"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PT4")),]$Patient <- "P4"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PN5" | orig.ident=="PT5")),]$Patient <- "P5"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PN6" | orig.ident=="PT6")),]$Patient <- "P6"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PT2N")),]$Patient <- "P7"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PT6N")),]$Patient <- "P8"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PT16P")),]$Patient <- "P9"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PT19P")),]$Patient <- "P10"

PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PN1" | orig.ident=="PN2" | orig.ident=="PN3" | orig.ident=="PN5" | orig.ident=="PN6")),]$Tissue <- "NATs"
PSCC@meta.data[rownames(subset(PSCC@meta.data,orig.ident=="PT1" | orig.ident=="PT2" | orig.ident=="PT3" | orig.ident=="PT4" | orig.ident=="PT5" | orig.ident=="PT6" | orig.ident=="PT6N" | orig.ident=="PT2N" | orig.ident=="PT16P" | orig.ident=="PT19P")),]$Tissue <- "Tumor"

PSCC@meta.data[rownames(subset(PSCC@meta.data,Patient=="P2" | Patient=="P3" | Patient=="P7" | Patient=="P8")),]$HPV <- "HPVN"
PSCC@meta.data[rownames(subset(PSCC@meta.data,Patient=="P1" | Patient=="P4" | Patient=="P5" | Patient=="P6" | Patient=="P9" | Patient=="P10")),]$HPV  <- "HPVP"

PSCC@meta.data[rownames(subset(PSCC@meta.data,Tissue=="NATs" & HPV=="HPVN")),]$Groups <- "NATs_HPVN"
PSCC@meta.data[rownames(subset(PSCC@meta.data,Tissue=="NATs" & HPV=="HPVP")),]$Groups <- "NATs_HPVP"
PSCC@meta.data[rownames(subset(PSCC@meta.data,Tissue=="Tumor" & HPV=="HPVN")),]$Groups <- "Tumor_HPVN"
PSCC@meta.data[rownames(subset(PSCC@meta.data,Tissue=="Tumor" & HPV=="HPVP")),]$Groups <- "Tumor_HPVP"

PSCC@meta.data[rownames(subset(PSCC@meta.data,Patient=="P1" | Patient=="P2" | Patient=="P3" | Patient=="P4" | Patient=="P5" | Patient=="P6")),]$Library  <- "3_scRNAseq"
PSCC@meta.data[rownames(subset(PSCC@meta.data,Patient=="P7" | Patient=="P8" | Patient=="P9" | Patient=="P10")),]$Library <- "5_scRNAseq"

VlnPlot(PSCC, features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.00,group.by = "orig.ident")
VlnPlot(subset(PSCC,Library=="3_scRNAseq"), features = c("nFeature_RNA","nCount_RNA","percent.mt"), ncol = 3,pt.size = 0.00,group.by = "orig.ident")

### ribosom and mitochondrion module scoring
ribosom.genes  <- rownames(PSCC)[grep("^RP[SL]",rownames(PSCC))] # 100 genes
mitochondrion.genes <-  rownames(PSCC)[grep("^MT-",rownames(PSCC))] # 13 genes
ribosome_mitochondrion <- c(ribosom.genes,mitochondrion.genes) #113 genes
C <-GetAssayData(object = PSCC, slot = "counts")
percent.mito <- Matrix::colSums(C[mitochondrion.genes,])/Matrix::colSums(C)*100
percent.ribo <- Matrix::colSums(C[ribosom.genes,])/Matrix::colSums(C)*100
PSCC <- AddMetaData(PSCC, percent.mito, col.name = "percent.mito")
PSCC <- AddMetaData(PSCC, percent.ribo, col.name = "percent.ribo")

### rds save
#saveRDS(PSCC.raw,"rds/PSCC.raw.rds") 25312 features, 148696 cells, title: orig.ident nCount_RNA nFeature_RNA  percent.mt
#saveRDS(PSCC,"rds/PSCC.rds") # 25312 features, 119351 cells, title: orig.ident nCount_RNA nFeature_RNA percent.mt Patient Tissue HPV Groups percent.mito percent.ribo Library
#saveRDS(PSCC[setdiff(rownames(PSCC),ribosome_mitochondrion),],"rds/PSCC.rmMtRb.rds") # 25312-113=25199 features, 119351 cells

### extract and analysis all cells from 10X3 (NATs and tumor) and 10X5 (tumor) scRNAseq data
rm(list=ls())
PSCC <- readRDS("rds/PSCC.rmMtRb.rds")

### define colors panel
colors.FeaturePlot <- c("lightgrey","red")
colors.umap <- paletteer_d("ggsci::default_igv")
colors.Patient <- paletteer_d("ggthemes::Tableau_10")
colors.Tissue <- c("#009966FF","#D60047FF")

###  Integrate data
batch.list <- SplitObject(PSCC,split.by = "Patient")
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures =2000)
})
features <- SelectIntegrationFeatures(object.list = batch.list)
anchors <- FindIntegrationAnchors(object.list = batch.list, anchor.features = features)
PSCC.combined <- IntegrateData(anchorset = anchors)

### UMAP analysis
DefaultAssay(PSCC.combined) <- "integrated"
PSCC.combined <- ScaleData(PSCC.combined,features = rownames(PSCC.combined))
PSCC.combined <- RunPCA(PSCC.combined,features = VariableFeatures(PSCC.combined))

PSCC.combined <- JackStraw(PSCC.combined,dims = 50)
PSCC.combined <- ScoreJackStraw(PSCC.combined , dims = 1:50)
JackStrawPlot(PSCC.combined, dims = 1:50)
ElbowPlot(PSCC.combined)

#PSCC.combined@meta.data <- PSCC.combined@meta.data[c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Patient","Tissue","HPV","Groups","Library")]
#PSCC.combined <- RunTSNE(PSCC.combined, dims = 1:10)
PSCC.combined <- RunUMAP(PSCC.combined, dims = 1:10)
PSCC.combined <- FindNeighbors(PSCC.combined,reduction = "pca",dims = 1:10)
PSCC.combined <- FindClusters(PSCC.combined, resolution =0.5)

### cell labels
Idents(PSCC.combined) <- "integrated_snn_res.0.5"
rank <- c(1,2,3,4,6,7,8,9,11,12,13,14,16,17,18,19,20,15,5,0,10)
levels(PSCC.combined) <- rank

new.cluster.ids <- c("C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21")
names(new.cluster.ids) <- levels(PSCC.combined)
PSCC.combined <- RenameIdents(PSCC.combined,new.cluster.ids)
PSCC.combined@meta.data["Cluster"] <- data.frame(Idents(PSCC.combined))[1]

new.cluster.ids <- c(rep("Immune cells",17),rep("Epithelial cells",1),rep("Endothelial cells",1),rep("CAFs",2))
names(new.cluster.ids) <- levels(PSCC.combined)
PSCC.combined <- RenameIdents(PSCC.combined,new.cluster.ids)
PSCC.combined@meta.data["CellTypes"] <- data.frame(Idents(PSCC.combined))[1]

### plot
Idents(PSCC.combined) <- "CellTypes"
DimPlot(PSCC.combined,reduction = "umap",raster = TALSE, label = T, repel = T,cols = colors.umap) # Idents(PSCC.combined) <- "CellTypes"
#DimPlot(PSCC.combined,reduction = "umap",raster = TALSE, label = T, repel = T, cols = colors.umap) # Idents(PSCC.combined) <- "Cluster"
#DimPlot(PSCC.combined,reduction = "umap",raster = TALSE, label = T, repel = T, cols = colors.Patient) # Idents(PSCC.combined) <- "Patient"
#DimPlot(PSCC.combined,reduction = "umap",raster = TALSE, label = T, cols = colors.Tissue) # Idents(PSCC.combined) <- "Library"

DefaultAssay(PSCC.combined) <- "RNA"
FeaturePlot(PSCC.combined,c("PTPRC","KRT5","PECAM1","COL6A2"),order = TRUE,raster = T,cols = colors.FeaturePlot)
FeaturePlot(PSCC.combined,c("ACTB","B2M","GAPDH"),order = TRUE,raster = T,cols = colors.FeaturePlot)

MySeuratWrappers::VlnPlot(PSCC.combined,features = c("PTPRC","KRT5","VWF","COL6A2"),pt.size =0,stack=T,x.lab = "",y.lab = "",direction = "horizontal",cols = colors.umap)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

### rds save
#saveRDS(PSCC.combined,"rds/PSCC.rmMtRb.combinedPatient.PC10.res05.rds")
