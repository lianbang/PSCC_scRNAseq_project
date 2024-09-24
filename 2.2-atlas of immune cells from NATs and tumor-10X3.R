
library(SeuratObject) # v4.1.3
library(Seurat) # v4.3.0
library(MySeuratWrappers)
library(ggplot2)
library(paletteer)

### Rscirpt for Figure 1I, 1J, 2A, 2B, S2F, S2G

### define colors panel
colors.FeaturePlot <- c("lightgrey","red")
colors.umap <- paletteer_d("ggsci::default_igv")
colors.Patient <- paletteer_d("ggthemes::Tableau_10")
colors.Tissue <- c("#009966FF","#D60047FF")
colors.HPV <-  c("#E6BB99","#C37A78")

### extract immune cells and reanalysis from NATs and tumor 10X3 scRNAseq data
rm(list=ls())
PSCC10X3 <- readRDS("rds/PSCC10X3.rmMtRb.combinedPatient.PC10.res01.rds")

PSCC10X3 <- subset(PSCC10X3, CellTypes=="Immune cells")
Immu.data <- GetAssayData(object = PSCC10X3, slot = "counts")

Immu10X3 <- CreateSeuratObject(counts = Immu.data, min.cells = 10, min.features = 200, project = "Immu10X3")
Immu10X3@meta.data <- PSCC@meta.data[rownames(Immu10X3@meta.data),]

###  Integrate data
batch.list <- SplitObject(Immu10X3,split.by = "Patient")
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures =2000)
})
features <- SelectIntegrationFeatures(object.list = batch.list)
anchors <- FindIntegrationAnchors(object.list = batch.list, anchor.features = features)
Immu10X3.combined <- IntegrateData(anchorset = anchors)

### UMAP analysis
DefaultAssay(Immu10X3.combined) <- "integrated"
Immu10X3.combined <- ScaleData(Immu10X3.combined,features = rownames(Immu10X3.combined))
Immu10X3.combined <- RunPCA(Immu10X3.combined,features = VariableFeatures(Immu10X3.combined))

Immu10X3.combined <- JackStraw(Immu10X3.combined,dims = 50)
Immu10X3.combined <- ScoreJackStraw(Immu10X3.combined , dims = 1:50)
JackStrawPlot(Immu10X3.combined, dims = 1:20)
ElbowPlot(Immu10X3.combined)

Immu10X3.combined@meta.data <- Immu10X3.combined@meta.data[c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Patient","Tissue","HPV","Groups","Library")]
#Immu10X3.combined <- RunTSNE(Immu10X3.combined, dims = 1:20)
Immu10X3.combined <- RunUMAP(Immu10X3.combined, dims = 1:20)
Immu10X3.combined <- FindNeighbors(Immu10X3.combined,reduction = "pca",dims = 1:20)
Immu10X3.combined <- FindClusters(Immu10X3.combined, resolution =0.5)

### cell labels
Idents(Immu10X3.combined) <- "integrated_snn_res.0.5"
rank <- c(5,6,3,2,11,1,8,16,19,10,4,0,7,13,9,17,18,14,12,15)
levels(Immu10X3.combined) <- rank

new.cluster.ids <- c("C1-B-MS4A1","C2-Plasma-SDC1", "C3-CD4T-FOXP3","C4-CD4T-CD40LG","C5-CD4T-CXCL13","C6-CD8T-GZMK","C7-CD8T-GNLY","C8-CD8T-MKI67","C9-PreT-PTCRA","C10-NK-KLRC1",
                     "C11-Mono-FCN1","C12-Macro-SELENOP","C13-Macro-SPP1","C14-Macro-C1QA","C15-DC-FCER1A","C16-DC-LAMP3","C17-DC-CD1A-MKI67","C18-Neutrophils-HCAR3","C19-Mast-MS4A2","C20-nonImmune-EPS8")
names(new.cluster.ids) <- levels(Immu10X3.combined)
Immu10X3.combined <- RenameIdents(Immu10X3.combined,new.cluster.ids)
Immu10X3.combined@meta.data["Anno"] <- data.frame(Idents(Immu10X3.combined))[1]

id <- c()
for (i in new.cluster.ids){
  id <- c(id, strsplit(i[1],split = '-')[[1]][1])
}
names(id) <- levels(Immu10X3.combined)
Immu10X3.combined <- RenameIdents(Immu10X3.combined,id)
Immu10X3.combined@meta.data["Cluster"] <- data.frame(Idents(Immu10X3.combined))[1]

new.cluster.ids <- c(rep("TIL-Bs",2),rep("TIL-Ts",7),rep("NK",1),rep("Myeloid",9),rep("nonImmune",1))
names(new.cluster.ids) <- levels(Immu10X3.combined)
Immu10X3.combined <- RenameIdents(Immu10X3.combined,new.cluster.ids)
Immu10X3.combined@meta.data["CellTypes"] <- data.frame(Idents(Immu10X3.combined))[1]

Idents(Immu10X3.combined) <- "Anno"
new.cluster.ids <- seq(1,20)
names(new.cluster.ids) <- levels(Immu10X3.combined)
Immu10X3.combined <- RenameIdents(Immu10X3.combined,new.cluster.ids)
Immu10X3.combined@meta.data["Rank"] <- data.frame(Idents(Immu10X3.combined))[1]

### Find Cluster Markers
Idents(Immu10X3.combined) <- "Cluster"
DefaultAssay(Immu10X3.combined) <- "RNA"
subtype.markers <- FindAllMarkers(Immu10X3.combined,logfc.threshold = 0.5,only.pos = TRUE) # Cluster idents
subset(subtype.markers,cluster=="C14")
subset(subtype.markers,gene=="HCAR3")

### plot
Idents(Immu10X3.combined) <- "Cluster"
DimPlot(Immu10X3.combined,reduction = "umap",order=T, raster = T, label = T,cols = colors.umap,repel = T,label.size = 6,)
DimPlot(Immu10X3.combined,reduction = "umap",order=T, raster = T, label = T,cols = colors.umap,repel = T,label.size = 6,split.by = "Tissue")
DimPlot(Immu10X3.combined,reduction = "umap",order=T, raster = T,label = T,cols = colors.umap,repel = T,group.by = "CellTypes",label.size = 4)
DimPlot(Immu10X3.combined,reduction = "umap",order=T, label = F,cols = colors.Patient,group.by = "Patient",pt.size = 0.1)
DimPlot(Immu10X3.combined,reduction = "umap",order=T, label = F,cols = colors.Tissue,group.by = "Tissue",pt.size = 0.1)

DefaultAssay(Immu10X3.combined) <- "RNA"
genes <- c("MS4A1","SDC1","CD3E","FOXP3","CD40LG","CXCL13","CD8A","GZMK","GNLY","MKI67","PTCRA","KLRC1","LYZ","CD86","FCN1","SELENOP","SPP1","C1QA","FCER1A","LAMP3","CD1A","HCAR3","MS4A2","EPS8")
MySeuratWrappers::VlnPlot(Immu10X3.combined,features = genes, pt.size =0,stack=T,x.lab = "",y.lab = "",direction = "horizontal",cols = colors.umap,line.size = 1)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

### relative proportion
Number.Patient_CellTypes <- xtabs(~orig.ident+CellTypes,data=Immu10X3.combined@meta.data)[,-5]
Percent.Patient_CellTypes <- data.frame(c(Number.Patient_CellTypes["PN1",],Number.Patient_CellTypes["PN2",],Number.Patient_CellTypes["PN3",],Number.Patient_CellTypes["PN5",],Number.Patient_CellTypes["PN6",],
                                       Number.Patient_CellTypes["PT1",],Number.Patient_CellTypes["PT2",],Number.Patient_CellTypes["PT3",],Number.Patient_CellTypes["PT4",],Number.Patient_CellTypes["PT5",],Number.Patient_CellTypes["PT6",]),# rownames
                                     as.character(rep(colnames(Number.Patient_CellTypes),dim(Number.Patient_CellTypes)[1])), #dim(Number.Patient_CellTypes)[1]
                                     c(rep("PN1",dim(Number.Patient_CellTypes)[2]),rep("PN2",dim(Number.Patient_CellTypes)[2]),rep("PN3",dim(Number.Patient_CellTypes)[2]),rep("PN5",dim(Number.Patient_CellTypes)[2]),rep("PN6",dim(Number.Patient_CellTypes)[2]),
                                       rep("PT1",dim(Number.Patient_CellTypes)[2]),rep("PT2",dim(Number.Patient_CellTypes)[2]),rep("PT3",dim(Number.Patient_CellTypes)[2]),rep("PT4",dim(Number.Patient_CellTypes)[2]),rep("PT5",dim(Number.Patient_CellTypes)[2]),rep("PT6",dim(Number.Patient_CellTypes)[2]))) #dim(Number.Patient_CellTypes)[2]
colnames(Percent.Patient_CellTypes) <- c("Number","CellTypes","Patient")

rank <- c("TIL-Bs","TIL-Ts","NK","Myeloid")
Percent.Patient_CellTypes$CellTypes = factor(Percent.Patient_CellTypes$CellTypes, levels =rank)

Percent.Patient_CellTypes.NATs <- subset(Percent.Patient_CellTypes,Patient=="PN1" | Patient=="PN2" | Patient=="PN3" | Patient=="PN5" | Patient=="PN6")
Percent.Patient_CellTypes.Tumor <- subset(Percent.Patient_CellTypes,Patient=="PT1" | Patient=="PT2" | Patient=="PT3" | Patient=="PT4" | Patient=="PT5" | Patient=="PT6")
Percent.Patient_CellTypes.Tumor$HPV_Patient <- "HPV_Patient"
Percent.Patient_CellTypes.Tumor$HPV_Patient <- c(rep("HPVP_P1",4),rep("HPVN_P2",4),rep("HPVN_P3",4),rep("HPVP_P4",4),rep("HPVP_P5",4),rep("HPVP_P6",4))

ggplot(data = Percent.Patient_CellTypes.Tumor, mapping =aes(x = HPV_Patient, y = Number,fill=CellTypes)) +
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

#Number.Tissue_CellTypes <- xtabs(~Tissue+CellTypes,data=Immu10X3.combined@meta.data)[,-5]
#Percent.Tissue_CellTypes <- data.frame(c(Number.Tissue_CellTypes["NATs",],Number.Tissue_CellTypes["Tumor",]),# rownames
#                                        as.character(rep(colnames(Number.Tissue_CellTypes),dim(Number.Tissue_CellTypes)[1])), #dim(Number.Tissue_CellTypes)[1]
#                                        c(rep("NATs",dim(Number.Tissue_CellTypes)[2]),rep("Tumor",dim(Number.Tissue_CellTypes)[2]))) #dim(Number.Tissue_CellTypes)[2]
#colnames(Percent.Tissue_CellTypes) <- c("Number","CellTypes","Tissue")

#rank <- c("TIL-Bs","TIL-Ts","NK","Myeloid")
#Percent.Tissue_CellTypes$CellTypes = factor(Percent.Tissue_CellTypes$CellTypes, levels =rank)

#ggplot(data = Percent.Tissue_CellTypes, mapping =aes(x = Tissue, y = Number,fill=CellTypes)) +
#  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
#  scale_fill_manual(values = colors.umap)+
#  ylab("Proportion")+
#geom_text(aes(label=number), position=position_stack(vjust=0.5)) +
#coord_flip()+ # 反转
#xlab("Tissue")+
#  theme_classic()+
#  theme(axis.title.x = element_blank(),
#        axis.title.y = element_text(colour="black", size=10),
#        axis.text.x = element_text(angle =0,size=10,colour="Black",hjust = 0.5,vjust=1),
#        axis.text.y = element_text(size=10,colour="Black"),
#        legend.title=element_text(size=10))# 图例名称 字体大小


### extract immune cells excluding for C20 from tumor tissue
Immu10X3.combined.tumor <- subset(Immu10X3.combined,Cluster!="C20")
Immu10X3.combined.tumor <- subset(Immu10X3.combined.tumor,Tissue=="Tumor")
DimPlot(Immu10X3.combined.tumor,reduction = "umap",order=T, raster = T, label = T,cols = colors.umap,repel = T,label.size = 6,)
DimPlot(Immu10X3.combined.tumor,reduction = "umap",order=T, raster = T, label = T,cols = colors.umap,repel = T,label.size = 6,split.by = "HPV")

### relative proportion
Number.HPV_Cluster <- xtabs(~HPV+Cluster,data=Immu10X3.combined.tumor@meta.data)
Percent.HPV_Cluster <- data.frame(c(Number.HPV_Cluster["HPVN",],Number.HPV_Cluster["HPVP",]),# rownames
                                  as.character(rep(colnames(Number.HPV_Cluster),dim(Number.HPV_Cluster)[1])), #dim(Number.HPV_Cluster)[1]
                                  c(rep("HPVN",dim(Number.HPV_Cluster)[2]),rep("HPVP",dim(Number.HPV_Cluster)[2]))) #dim(Number.HPV_Cluster)[2]

colnames(Percent.HPV_Cluster) <- c("Number","Cluster","HPV")
id <- colnames(Number.HPV_Cluster)
Percent.HPV_Cluster$Cluster = factor(Percent.HPV_Cluster$Cluster, levels =id) # id, C1-C19

ggplot(data = Percent.HPV_Cluster, mapping =aes(x = Cluster, y = Number,fill= HPV)) +
  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
  scale_fill_manual(values = colors.HPV)+
  xlab("Clusters")+
  ylab("Relative proportion")+
  #coord_flip()+ # 反转
  theme_classic(12)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black", size=12),
        axis.text.x = element_text(angle =45,size=12,colour="Black",hjust = 1,vjust=1),
        axis.text.y = element_text(size=12,colour="Black"),
        legend.title=element_text(size=12))# 图例名称 字体大小

### rds save
#saveRDS(Immu10X3.combined,"rds/Immu10X3.rmMtRb.combinedPatient.PC20.res05.rds")
#saveRDS(subtype.markers,"rds/subtype.markers.Immu10X3.rmMtRb.combinedPatient.PC20.res05.rds")
