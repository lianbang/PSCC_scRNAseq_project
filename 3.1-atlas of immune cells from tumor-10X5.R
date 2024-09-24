
library(SeuratObject) # v4.1.3
library(Seurat) # v4.3.0
library(MySeuratWrappers)
library(ggplot2)
library(paletteer)

### define colors panel
colors.FeaturePlot <- c("lightgrey","red")
colors.umap <- paletteer_d("ggsci::default_igv")
colors.Patient <- paletteer_d("ggthemes::Tableau_10")
colors.Tissue <- c("#009966FF","#D60047FF")
colors.HPV <-  c("#E6BB99","#C37A78")

### extract and analysis sorted immune cells from tumor 10X3 scRNAseq data
PSCC <- readRDS("rds/PSCC.rmMtRb.rds")
#PSCC <- readRDS("rds/PSCC.rmMtRb.combinedPatient.PC10.res05.rds")
PSCC <- subset(PSCC, Library=="5_scRNAseq")

PSCC.data <-GetAssayData(object = PSCC, slot = "counts")
PSCC10X5 <- CreateSeuratObject(counts = PSCC.data, min.cells = 10, min.features = 200, project = "10X5")
PSCC10X5@meta.data <- PSCC@meta.data[rownames(PSCC10X5@meta.data),]

###  Integrate data
batch.list <- SplitObject(PSCC10X5,split.by = "Patient")
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures =2000)
})
features <- SelectIntegrationFeatures(object.list = batch.list)
anchors <- FindIntegrationAnchors(object.list = batch.list, anchor.features = features)
PSCC10X5.combined <- IntegrateData(anchorset = anchors)

### UMAP analysis
DefaultAssay(PSCC10X5.combined) <- "integrated"
PSCC10X5.combined <- ScaleData(PSCC10X5.combined,features = rownames(PSCC10X5.combined))
PSCC10X5.combined <- RunPCA(PSCC10X5.combined,features = VariableFeatures(PSCC10X5.combined))

PSCC10X5.combined <- JackStraw(PSCC10X5.combined,dims = 50)
PSCC10X5.combined <- ScoreJackStraw(PSCC10X5.combined , dims = 1:50)
JackStrawPlot(PSCC10X5.combined, dims = 1:50)
ElbowPlot(PSCC10X5.combined)

#PSCC10X5.combined@meta.data <- PSCC10X5.combined@meta.data[c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Patient","Tissue","HPV","Groups","Library")]
#PSCC10X5.combined <- RunTSNE(PSCC10X5.combined, dims = 1:10)
PSCC10X5.combined <- RunUMAP(PSCC10X5.combined, dims = 1:10)
PSCC10X5.combined <- FindNeighbors(PSCC10X5.combined,reduction = "pca",dims = 1:10)
PSCC10X5.combined <- FindClusters(PSCC10X5.combined, resolution =0.1)

### cell labels
Idents(PSCC10X5.combined) <- "integrated_snn_res.0.1"
rank <- c(0,1,2,3,5,4,7,6)
levels(PSCC10X5.combined) <- rank

new.cluster.ids <- c(rep("TIL-Ts",5),rep("TIL-Bs",2),"Myeloid")
names(new.cluster.ids) <- levels(PSCC10X5.combined)
PSCC10X5.combined <- RenameIdents(PSCC10X5.combined,new.cluster.ids)
PSCC10X5.combined@meta.data["CellTypes"] <- data.frame(Idents(PSCC10X5.combined))[1]

Idents(PSCC10X5.combined) <- "integrated_snn_res.0.1"
levels(PSCC10X5.combined) <- rank
names(new.cluster.ids) <- levels(PSCC10X5.combined)
PSCC10X5.combined <- RenameIdents(PSCC10X5.combined,new.cluster.ids)
PSCC10X5.combined@meta.data["Cluster"] <- data.frame(Idents(PSCC10X5.combined))[1]

### Find Cluster Markers
Idents(PSCC10X5.combined) <- "Cluster"
DefaultAssay(PSCC10X5.combined) <- "RNA"
subtype.markers <- FindAllMarkers(PSCC10X5.combined,logfc.threshold = 0.5,only.pos = TRUE) # Cluster idents
subset(subtype.markers,cluster=="C1")
subset(subtype.markers,gene=="KLRC1")

### plot
Idents(PSCC10X5.combined) <- "Cluster"
DimPlot(PSCC10X5.combined,reduction = "umap",order=T, raster = T,label = T,cols = colors.umap,repel = T)

DefaultAssay(PSCC10X5.combined) <- "RNA"
FeaturePlot(PSCC10X5.combined,c("PTPRC","CD3E","MS4A1","SDC1","LYZ"),order = F,cols = colors.FeaturePlot, raster = T, ncol = 3)
DotPlot(PSCC10X5.combined,features = c("PTPRC","CD3E","MS4A1","SDC1","LYZ"))

MySeuratWrappers::VlnPlot(PSCC10X5.combined,features = c("PTPRC","CD3E","MS4A1","SDC1","LYZ"),pt.size =0,stack=T,x.lab = "",y.lab = "",direction = "horizontal")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


### relative proportion
Number.Patient_CellTypes <- xtabs(~orig.ident+CellTypes,data=PSCC10X5.combined@meta.data)
Percent.Patient_CellTypes <- data.frame(c(Number.Patient_CellTypes["PT16P",],Number.Patient_CellTypes["PT19P",],Number.Patient_CellTypes["PT2N",],Number.Patient_CellTypes["PT6N",]),# rownames
                                        as.character(rep(colnames(Number.Patient_CellTypes),dim(Number.Patient_CellTypes)[1])), #dim(Number.Patient_CellTypes)[1]
                                        c(rep("PT16P",dim(Number.Patient_CellTypes)[2]),rep("PT19P",dim(Number.Patient_CellTypes)[2]),rep("PT2N",dim(Number.Patient_CellTypes)[2]),rep("PT6N",dim(Number.Patient_CellTypes)[2]))) #dim(Number.Patient_CellTypes)[2]

colnames(Percent.Patient_CellTypes) <- c("Number","CellTypes","Patient")
rank <- c("TIL-Bs","TIL-Ts","Myeloid")
Percent.Patient_CellTypes$CellTypes = factor(Percent.Patient_CellTypes$CellTypes, levels =rank) #
Percent.Patient_CellTypes$HPV_Patient <- "HPV_Patient"
Percent.Patient_CellTypes$HPV_Patient <- c("HPVP_P9","HPVP_10","HPVN_P7","HPVN_P8")

ggplot(data = Percent.Patient_CellTypes, mapping =aes(x = HPV_Patient, y = Number,fill=CellTypes)) +
  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
  scale_fill_manual(values = colors.umap[c(1,2,4)])+
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

### rds save
#saveRDS(PSCC10X5.combined,"rds/PSCC10X5.rmMtRb.combinedPatient.PC10.res01.rds")
