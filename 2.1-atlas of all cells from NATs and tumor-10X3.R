
library(SeuratObject) # v4.1.3
library(Seurat) # v4.3.0
library(MySeuratWrappers)
library(ggplot2)
library(paletteer)

### Rscirpt for Figure 1B, 1C, S2D

### define colors panel
colors.FeaturePlot <- c("lightgrey","red")
colors.umap <- paletteer_d("ggsci::default_igv")
colors.Patient <- paletteer_d("ggthemes::Tableau_10")
colors.Tissue <- c("#009966FF","#D60047FF")

### extract and analysis all cells from NATs and tumor 10X3 scRNAseq data
rm(list=ls())
PSCC <- readRDS("rds/PSCC.rmMtRb.rds")
#PSCC <- readRDS("rds/PSCC.rmMtRb.combinedPatient.PC10.res05.rds")

PSCC <- subset(PSCC, Library=="3_scRNAseq")
PSCC.data <-GetAssayData(object = PSCC, slot = "counts")

PSCC10X3 <- CreateSeuratObject(counts = PSCC.data, min.cells = 10, min.features = 200, project = "10X3")
PSCC10X3@meta.data <- PSCC@meta.data[rownames(PSCC10X3@meta.data),]

###  Integrate data
batch.list <- SplitObject(PSCC10X3,split.by = "Patient")
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures =2000)
})
features <- SelectIntegrationFeatures(object.list = batch.list)
anchors <- FindIntegrationAnchors(object.list = batch.list, anchor.features = features)
PSCC10X3.combined <- IntegrateData(anchorset = anchors)

### UMAP analysis
DefaultAssay(PSCC10X3.combined) <- "integrated"
PSCC10X3.combined <- ScaleData(PSCC10X3.combined,features = rownames(PSCC10X3.combined))
PSCC10X3.combined <- RunPCA(PSCC10X3.combined,features = VariableFeatures(PSCC10X3.combined))

PSCC10X3.combined <- JackStraw(PSCC10X3.combined,dims = 50)
PSCC10X3.combined <- ScoreJackStraw(PSCC10X3.combined , dims = 1:50)
JackStrawPlot(PSCC10X3.combined, dims = 1:50)
ElbowPlot(PSCC10X3.combined)

#PSCC10X3.combined@meta.data <- PSCC10X3.combined@meta.data[c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Patient","Tissue","HPV","Groups","Library")]
#PSCC10X3.combined <- RunTSNE(PSCC10X3.combined, dims = 1:10)
PSCC10X3.combined <- RunUMAP(PSCC10X3.combined, dims = 1:10)
PSCC10X3.combined <- FindNeighbors(PSCC10X3.combined,reduction = "pca",dims = 1:10)
PSCC10X3.combined <- FindClusters(PSCC10X3.combined, resolution =0.1)

### cell labels
Idents(PSCC10X3.combined) <- "integrated_snn_res.0.1"
rank <- c(6,7,2,3,8,5,1,9,0,4)
levels(PSCC10X3.combined) <- rank

new.cluster.ids <- c("C1-B cells","C2-Plamsa cells","C3-T cells","C4-Myeloid cells","C5-Mast cells","C6-Epithelial cells","C7-Endothelial cells-1","C8-Endothelial cells-2","C9-iCAFs","C10-MyoCAFs")
names(new.cluster.ids) <- levels(PSCC10X3.combined)
PSCC10X3.combined <- RenameIdents(PSCC10X3.combined,new.cluster.ids)
PSCC10X3.combined@meta.data["Anno"] <- data.frame(Idents(PSCC10X3.combined))[1]

Idents(PSCC10X3.combined) <- "Anno"
id <- c()
for (i in new.cluster.ids){
  id <- c(id, strsplit(i[1],split = '-')[[1]][1])
}
names(id) <- levels(PSCC10X3.combined)
PSCC10X3.combined<- RenameIdents(PSCC10X3.combined,id)
PSCC10X3.combined@meta.data["Cluster"] <- data.frame(Idents(PSCC10X3.combined))[1]

new.cluster.ids <- c("Immune cells","Immune cells","Immune cells","Immune cells","Immune cells","Epithelial cells","Endothelial cells","Endothelial cells","CAFs","CAFs")
names(new.cluster.ids) <- levels(PSCC10X3.combined)
PSCC10X3.combined <- RenameIdents(PSCC10X3.combined,new.cluster.ids)
PSCC10X3.combined@meta.data["CellTypes"] <- data.frame(Idents(PSCC10X3.combined))[1]

Idents(PSCC10X3.combined) <- "Anno"
new.cluster.ids <- seq(1,10)
names(new.cluster.ids) <- levels(PSCC10X3.combined)
PSCC10X3.combined <- RenameIdents(PSCC10X3.combined,new.cluster.ids)
PSCC10X3.combined@meta.data["Rank"] <- data.frame(Idents(PSCC10X3.combined))[1]

### Find Cluster Markers
Idents(PSCC10X3.combined) <- "Cluster"
DefaultAssay(PSCC10X3.combined) <- "RNA"
subtype.markers <- FindAllMarkers(PSCC10X3.combined,logfc.threshold = 0.5,only.pos = TRUE) # Cluster idents
subset(subtype.markers,cluster=="C1")
subset(subtype.markers,gene=="ECSCR")

### Find CellTypes markers, including Immune cells, Epithelial cells, Endothelial cells, CAFs
Idents(PSCC10X3.combined) <- "CellTypes"
DefaultAssay(PSCC10X3.combined) <- "RNA"
CellTypes.markers <- FindAllMarkers(PSCC10X3.combined,logfc.threshold = 0.5,only.pos = TRUE) # CellTypes idents
subset(CellTypes.markers,cluster=="CAFs")
subset(CellTypes.markers,gene=="PTPRC")

### module scoring
CAFs.genes <-  c("LUM","COL1A1","COL3A1","DCN","CFD")
EC.genes <- c("ACKR1","RAMP2","SELE","VWF","PECAM1")

C <-GetAssayData(object = PSCC10X3.combined, slot = "counts")
CAFs.genes <- Matrix::colSums(C[CAFs.genes,])/Matrix::colSums(C)*100
EC.genes <- Matrix::colSums(C[EC.genes,])/Matrix::colSums(C)*100
PSCC10X3.combined <- AddMetaData(PSCC10X3.combined, CAFs.genes, col.name = "CAFs.genes")
PSCC10X3.combined <- AddMetaData(PSCC10X3.combined, EC.genes, col.name = "EC.genes")

### plot
Idents(PSCC10X3.combined) <- "CellTypes"
DimPlot(PSCC10X3.combined,reduction = "umap",order=T, label = T,cols = colors.umap,repel = T)
#DimPlot(PSCC10X3.combined,reduction = "umap",order=T, raster = T,label = T,cols = colors.umap,repel = T) # Idents(PSCC10X3.combined) <- "Cluster"
DimPlot(PSCC10X3.combined,reduction = "umap",order=T, label = F, cols = colors.Patient,group.by = "Patient",pt.size = 0.1)
DimPlot(PSCC10X3.combined,reduction = "umap",order=T, label = F, cols = colors.Tissue,group.by = "Tissue",pt.size = 0.1)

DefaultAssay(PSCC10X3.combined) <- "RNA"
FeaturePlot(PSCC10X3.combined,c("PTPRC","KRT5","PECAM1","COL6A2"),order = T,cols = colors.FeaturePlot, raster = T, ncol = 4,pt.size = 1)

DotPlot(PSCC10X3.combined,features = c("PTPRC","KRT5","PECAM1","RAMP2","SELE","VWF","ACKR1","COL6A2","COL1A1","COL3A1","LUM","CFD","DCN"),cols = colors.FeaturePlot) +coord_flip()+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y.left = element_blank(),
        axis.title.y = element_text(colour="black", size=12),
        axis.text.x = element_text(angle =45,size=12,colour="Black",hjust = 1,vjust=1),
        axis.text.y = element_text(angle =0,size=12,colour="Black",hjust = 0.5,vjust=1))

#DotPlot(PSCC10X3.combined,features = c("PTPRC","KRT5","PECAM1","COL6A2"))

#MySeuratWrappers::VlnPlot(PSCC10X3.combined,features = c("PTPRC","KRT5","CAFs.genes","EC.genes"),pt.size =0,stack=T,x.lab = "",y.lab = "",direction = "horizontal")+
#  theme(axis.text.x = element_blank(),
#        axis.ticks.x = element_blank())

top <- CellTypes.markers  %>% group_by(cluster) %>% top_n(20, avg_log2FC)
p <- DoHeatmap(object = PSCC10X3.combined,features = top$gene, label = TRUE, group.colors = colors.umap, assay = "integrated")+NoLegend() # "integrated" Assay

set.seed(42)
subobj <- subset(PSCC10X3.combined, downsample = 1000)
p <- DoHeatmap(object = subobj,features = top$gene, label = TRUE, group.colors = colors.umap, assay = "integrated")+NoLegend() # "integrated" Assay


### relative proportion
id <- c("Immune cells", "Epithelial cells", "Endothelial cells", "CAFs")
Number.Patient_Tissue <- xtabs(~orig.ident+CellTypes,data=PSCC10X3.combined@meta.data)
Percent.Patient_Tissue <- data.frame(c(Number.Patient_Tissue["PN1",],Number.Patient_Tissue["PN2",],Number.Patient_Tissue["PN3",],Number.Patient_Tissue["PN5",],Number.Patient_Tissue["PN6",],
                                       Number.Patient_Tissue["PT1",],Number.Patient_Tissue["PT2",],Number.Patient_Tissue["PT3",],Number.Patient_Tissue["PT4",],Number.Patient_Tissue["PT5",],Number.Patient_Tissue["PT6",]),# rownames
                                     as.character(rep(colnames(Number.Patient_Tissue),dim(Number.Patient_Tissue)[1])), #dim(Number.Patient_Tissue)[1]
                                     c(rep("PN1",dim(Number.Patient_Tissue)[2]),rep("PN2",dim(Number.Patient_Tissue)[2]),rep("PN3",dim(Number.Patient_Tissue)[2]),rep("PN5",dim(Number.Patient_Tissue)[2]),rep("PN6",dim(Number.Patient_Tissue)[2]),
                                       rep("PT1",dim(Number.Patient_Tissue)[2]),rep("PT2",dim(Number.Patient_Tissue)[2]),rep("PT3",dim(Number.Patient_Tissue)[2]),rep("PT4",dim(Number.Patient_Tissue)[2]),rep("PT5",dim(Number.Patient_Tissue)[2]),rep("PT6",dim(Number.Patient_Tissue)[2]))) #dim(Number.Patient_Tissue)[2]
colnames(Percent.Patient_Tissue) <- c("Number","CellTypes","Patient_Tissue")

Percent.Patient_Tissue$CellTypes = factor(Percent.Patient_Tissue$CellTypes, levels =id)

Percent.Patient_Tissue.NATs <- subset(Percent.Patient_Tissue,Patient_Tissue=="PN1" | Patient_Tissue=="PN2" | Patient_Tissue=="PN3" | Patient_Tissue=="PN5" | Patient_Tissue=="PN6")
Percent.Patient_Tissue.Tumor <- subset(Percent.Patient_Tissue,Patient_Tissue=="PT1" | Patient_Tissue=="PT2" | Patient_Tissue=="PT3" | Patient_Tissue=="PT4" | Patient_Tissue=="PT5" | Patient_Tissue=="PT6")

ggplot(data = Percent.Patient_Tissue.Tumor, mapping =aes(x = Patient_Tissue, y = Number,fill=CellTypes)) +
  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
  scale_fill_manual(values = colors.umap)+
  ylab("Proportion")+
  #geom_text(aes(label=number), position=position_stack(vjust=0.5)) +
  coord_flip()+ # 反转
  #xlab("Patient")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black", size=10),
        axis.text.x = element_text(angle =45,size=10,colour="Black",hjust = 1,vjust=1),
        axis.text.y = element_text(size=10,colour="Black"),
        legend.title=element_text(size=10))# 图例名称 字体大小


#Number.Tissue <- xtabs(~Tissue+CellTypes,data=PSCC10X3.combined@meta.data)
#Percent.Tissue <- data.frame(c(Number.Tissue["NATs",],Number.Tissue["Tumor",]),# rownames
#                                     as.character(rep(colnames(Number.Tissue),dim(Number.Tissue)[1])), #dim(Number.Tissue)[1]
#                                     c(rep("NATs",dim(Number.Tissue)[2]),rep("Tumor",dim(Number.Tissue)[2]))) #dim(Number.Tissue)[2]

#colnames(Percent.Tissue) <- c("Number","CellTypes","Tissue")
#Percent.Tissue$CellTypes = factor(Percent.Tissue$CellTypes, levels =id)

#ggplot(data = Percent.Tissue, mapping =aes(x = Tissue, y = Number, fill=CellTypes)) +
#  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
#  scale_fill_manual(values = colors.umap)+
#  ylab("Proportion")+
  #geom_text(aes(label=number), position=position_stack(vjust=0.5)) +
  #coord_flip()+ # 反转
  #xlab("Patient")+
#  theme_classic()+
#  theme(axis.title.x = element_blank(),
#        axis.title.y = element_text(colour="black", size=10),
#        axis.text.x = element_text(angle =0,size=10,colour="Black",hjust = 0.5,vjust=1),
#        axis.text.y = element_text(size=10,colour="Black"),
#        legend.title=element_text(size=10), # 图例名称 字体大小
#        #legend.margin=unit(10,"cm"),
#        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

### rds save
#saveRDS(PSCC10X3.combined,"rds/PSCC10X3.rmMtRb.combinedPatient.PC10.res01.rds") # 25349 features, 44811 cells
#saveRDS(subtype.markers,"rds/subtype.markers.PSCC10X3.rmMtRb.combinedPatient.PC10.res01.rds")
#saveRDS(CellTypes.markers,"rds/CellTypes.markers.PSCC10X3.rmMtRb.combinedPatient.PC10.res01.rds")
