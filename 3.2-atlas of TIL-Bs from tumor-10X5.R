
library(SeuratObject) # v4.1.3
library(Seurat) # v4.3.0
library(MySeuratWrappers)
library(ggplot2)
library(ggrepel)
library(paletteer)
library(scRepertoire)
library(ggsci)
library(clusterProfiler) # bitr
library(org.Hs.eg.db)

### Rscirpt for Figure 3A, 3B, S4A, S4B, S4C, S4G, S4E

### define colors panel
colors.FeaturePlot <- c("lightgrey","red")
colors.umap <- paletteer_d("ggsci::default_igv")
colors.Patient <- paletteer_d("ggthemes::Tableau_10")
colors.Tissue <- c("#009966FF","#D60047FF")
#colors.HPV <-  c("#E6BB99","#C37A78")
colors.HPV <-  c("#660099FF","#FF1463FF")
colors.BCR <- c("Yellow","#FF1463FF")

### Integrate single cell BCR contig
contig_BCR_P7 <- read.csv("cellranger-6.0.1/scBCR-seq/P7-BCR-filtered_contig_annotations.csv")
contig_BCR_P8 <- read.csv("cellranger-6.0.1/scBCR-seq/P8-BCR-filtered_contig_annotations.csv")
contig_BCR_P9 <- read.csv("cellranger-6.0.1/scBCR-seq/P9-BCR-filtered_contig_annotations.csv")
contig_BCR_P10 <- read.csv("cellranger-6.0.1/scBCR-seq/P10-BCR-filtered_contig_annotations.csv")
contig_BCR <- list(contig_BCR_P7,contig_BCR_P8,contig_BCR_P9,contig_BCR_P10)
BCR <- combineBCR(contig_BCR,samples = c("HPVN1","HPVN2","HPVP1","HPVP2"))
BCR.2 <- combineBCR(contig_BCR,samples = c("PT2N","PT6N","PT16P","PT19P"))

head(BCR$HPVN1) #head(BCR$HPVN1)
subsetContig(BCR, name = "sample", variables = c("HPVN1","HPVN2","HPVP1","HPVP2"))

BCR_P7_barcode <- unlist(lapply(BCR$HPVN1$barcode,function(x){strsplit(x,"_")[[1]][2]}))
BCR_P7_barcode <- paste("PT2N",BCR_P7_barcode,sep = "_")
BCR_P8_barcode <- unlist(lapply(BCR$HPVN2$barcode,function(x){strsplit(x,"_")[[1]][2]}))
BCR_P8_barcode <- paste("PT6N",BCR_P8_barcode,sep = "_")
BCR_P9_barcode <- unlist(lapply(BCR$HPVP1$barcode,function(x){strsplit(x,"_")[[1]][2]}))
BCR_P9_barcode <- paste("PT16P",BCR_P9_barcode,sep = "_")
BCR_P10_barcode <- unlist(lapply(BCR$HPVP2$barcode,function(x){strsplit(x,"_")[[1]][2]}))
BCR_P10_barcode <- paste("PT19P",BCR_P10_barcode,sep = "_")

### files and rds save
#write.table(BCR$HPVN1,"rds/BCR.HPVN1.P7.txt",sep = "\t",quote = FALSE) # 2395 cells
#write.table(BCR$HPVN2,"rds/BCR.HPVN2.P8.txt",sep = "\t",quote = FALSE) # 976 cells
#write.table(BCR$HPVP1,"rds/BCR.HPVP1.P9.txt",sep = "\t",quote = FALSE) # 3312 cells
#write.table(BCR$HPVP2,"rds/BCR.HPVP2.P10.txt",sep = "\t",quote = FALSE) # 3109 cells
#saveRDS(BCR,"rds/BCR.scRepertoire.rds")
#saveRDS(BCR.2,"rds/BCR.2.scRepertoire.rds")

### extract and analysis TIL-Bs paired with BCR tumor 10X3 scRNAseq data
PSCC10X5.combined <- readRDS("rds/PSCC10X5.rmMtRb.combinedPatient.PC10.res01.rds")

PSCC10X5.BCR <- PSCC10X5.combined[,c(BCR_P7_barcode,BCR_P8_barcode,BCR_P9_barcode,BCR_P10_barcode)]
BCR10X5.data <-GetAssayData(object = PSCC10X5.BCR, slot = "counts")

BCR10X5 <- CreateSeuratObject(counts = BCR10X5.data, min.cells = 10, min.features = 200, project = "BCR10X5")
BCR10X5@meta.data <- PSCC10X5.BCR@meta.data[rownames(BCR10X5@meta.data),]
BCR10X5@meta.data <- BCR10X5@meta.data[c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Patient","Tissue","HPV","Groups","Library")]

###  Integrate data
batch.list <- SplitObject(BCR10X5,split.by = "Patient")
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures =2000)
})
features <- SelectIntegrationFeatures(object.list = batch.list)
anchors <- FindIntegrationAnchors(object.list = batch.list, anchor.features = features)
BCR10X5.combined <- IntegrateData(anchorset = anchors)

### UMAP analysis
DefaultAssay(BCR10X5.combined) <- "integrated"
BCR10X5.combined <- ScaleData(BCR10X5.combined,features = rownames(BCR10X5.combined))
BCR10X5.combined <- RunPCA(BCR10X5.combined,features = VariableFeatures(BCR10X5.combined))

BCR10X5.combined <- JackStraw(BCR10X5.combined,dims = 50)
BCR10X5.combined <- ScoreJackStraw(BCR10X5.combined , dims = 1:50)
JackStrawPlot(BCR10X5.combined, dims = 1:20)
ElbowPlot(BCR10X5.combined)

#BCR10X5.combined <- RunTSNE(BCR10X5.combined, dims = 1:10)
BCR10X5.combined <- RunUMAP(BCR10X5.combined, dims = 1:10)
BCR10X5.combined <- FindNeighbors(BCR10X5.combined,reduction = "pca",dims = 1:10)
BCR10X5.combined <- FindClusters(BCR10X5.combined, resolution =0.1)

### plot
DimPlot(BCR10X5.combined,reduction = "umap",order=T, raster = T, label = T,cols = colors.umap,repel = T,label.size = 6,)

DefaultAssay(BCR10X5.combined) <- "RNA"
genes <- c("PTPRC","CD19","MS4A1","SDC1","CD38","CD3E")
MySeuratWrappers::VlnPlot(BCR10X5.combined,features = genes, pt.size =0,stack=T,x.lab = "",y.lab = "",direction = "horizontal",cols = colors.umap,line.size = 1)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

FeaturePlot(BCR10X5.combined,genes,order = F,cols = colors.FeaturePlot, raster = F, ncol = 3)
DotPlot(BCR10X5.combined,features = genes)

### rds save
#saveRDS(BCR10X5.combined,"rds/BCR10X5.combinedPatient.PC20.res01.allBCR.8102cells.rds")

### extract and analysis TIL-Bs paired with BCR tumor 10X3 scRNAseq data (exclude for T cells)
DefaultAssay(BCR10X5.combined) <- "RNA"
BCR10X5.combined.C04 <- subset(BCR10X5.combined, integrated_snn_res.0.1=="0" | integrated_snn_res.0.1=="4")
BCR10X5.combined.C04 <- subset(BCR10X5.combined.C04, CD3D==0 & CD3E==0 & CD3G==0 & CD4==0 & CD8A==0 & CD8B==0)
TILBs.BCR10X5.data <- GetAssayData(object = BCR10X5.combined.C04, slot = "counts")

TILBs.BCR10X5 <- CreateSeuratObject(counts = TILBs.BCR10X5.data, min.cells = 10, min.features = 200, project = "TILBs.BCR10X5")
TILBs.BCR10X5@meta.data <- BCR10X5.combined.C04@meta.data[rownames(TILBs.BCR10X5@meta.data),]
TILBs.BCR10X5@meta.data <- TILBs.BCR10X5@meta.data[c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Patient","Tissue","HPV","Groups","Library")]

###  Integrate data
batch.list <- SplitObject(TILBs.BCR10X5,split.by = "HPV")
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures =2000)
})
features <- SelectIntegrationFeatures(object.list = batch.list)
anchors <- FindIntegrationAnchors(object.list = batch.list, anchor.features = features)
TILBs.BCR10X5.combined <- IntegrateData(anchorset = anchors)

### UMAP analysis
DefaultAssay(TILBs.BCR10X5.combined) <- "integrated"
TILBs.BCR10X5.combined <- ScaleData(TILBs.BCR10X5.combined,features = rownames(TILBs.BCR10X5.combined))
TILBs.BCR10X5.combined <- RunPCA(TILBs.BCR10X5.combined,features = VariableFeatures(TILBs.BCR10X5.combined))

TILBs.BCR10X5.combined <- JackStraw(TILBs.BCR10X5.combined,dims = 50)
TILBs.BCR10X5.combined <- ScoreJackStraw(TILBs.BCR10X5.combined , dims = 1:50)
JackStrawPlot(TILBs.BCR10X5.combined, dims = 1:20)
ElbowPlot(TILBs.BCR10X5.combined)

#TILBs.BCR10X5.combined <- RunTSNE(TILBs.BCR10X5.combined, dims = 1:6)
TILBs.BCR10X5.combined <- RunUMAP(TILBs.BCR10X5.combined, dims = 1:6)
TILBs.BCR10X5.combined <- FindNeighbors(TILBs.BCR10X5.combined,reduction = "pca",dims = 1:6)
TILBs.BCR10X5.combined <- FindClusters(TILBs.BCR10X5.combined, resolution =0.4)

### cell labels
Idents(TILBs.BCR10X5.combined) <- "integrated_snn_res.0.4"
rank <- c(4,0,2,1,5,6,3)
levels(TILBs.BCR10X5.combined) <- rank

new.cluster.ids <- c("C1-Naive B","C2-Memory B(EGR1-)","C2-Memory B(EGR1-)", "C3-Memory B(EGR1+)","C4-GC B","C5-Plasmablast","C6-Plasma cells")
names(new.cluster.ids) <- levels(TILBs.BCR10X5.combined)
TILBs.BCR10X5.combined <- RenameIdents(TILBs.BCR10X5.combined,new.cluster.ids)
TILBs.BCR10X5.combined@meta.data["Anno"] <- data.frame(Idents(TILBs.BCR10X5.combined))[1]

Idents(TILBs.BCR10X5.combined) <- "Anno"
new.cluster.ids <- c("C1","C2","C3","C4","C5","C6")
names(new.cluster.ids) <- levels(TILBs.BCR10X5.combined)
TILBs.BCR10X5.combined <- RenameIdents(TILBs.BCR10X5.combined,new.cluster.ids)
TILBs.BCR10X5.combined@meta.data["Cluster"] <- data.frame(Idents(TILBs.BCR10X5.combined))[1]

Idents(TILBs.BCR10X5.combined) <- "Anno"
new.cluster.ids <- seq(1,6)
names(new.cluster.ids) <- levels(TILBs.BCR10X5.combined)
TILBs.BCR10X5.combined <- RenameIdents(TILBs.BCR10X5.combined,new.cluster.ids)
TILBs.BCR10X5.combined@meta.data["Rank"] <- data.frame(Idents(TILBs.BCR10X5.combined))[1]

### Find Cluster Markers
Idents(TILBs.BCR10X5.combined) <- "Cluster"
DefaultAssay(TILBs.BCR10X5.combined) <- "RNA"
subtype.markers <- FindAllMarkers(TILBs.BCR10X5.combined,logfc.threshold = 0.5,only.pos = TRUE)
subset(subtype.markers,gene=="EGR1")

### plot
DimPlot(TILBs.BCR10X5.combined,reduction = "umap",order=T, raster = T, label = T,cols = colors.umap,repel = T,label.size = 6,pt.size = 0.2) # Idents(TILBs.BCR10X5.combined) <- "Cluster"

genes <- c("CD19","MS4A1","IGHD","IGHM","CD24","CD27","EGR1","CD38","SDC1","MKI67","BCL6","AICDA","MME","SUGCT","IRF4")
MySeuratWrappers::VlnPlot(TILBs.BCR10X5.combined,features = genes, pt.size =0,stack=T,x.lab = "",y.lab = "",direction = "horizontal",cols = colors.umap,line.size = 1)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

genes <- c("CD19","MS4A1","IGHD","IGHM","CD24","CD27","EGR1","CD38","SDC1","MKI67","BCL6","AICDA","MME","SUGCT")
MySeuratWrappers::VlnPlot(TILBs.BCR10X5.combined,features = genes, pt.size =0,stack=T,x.lab = "",y.lab = "",direction = "horizontal",cols = colors.umap,line.size = 1)+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#FeaturePlot(TILBs.BCR10X5.combined,features = c("PTPRC","CD19","MS4A1","AICDA"), order = T ,cols = colors.FeaturePlot, pt.size = 0.2)
#FeaturePlot(TILBs.BCR10X5.combined,features = c("PTPRC","CD19","MS4A1","CD79A","CD79B","IGHD","IGHM","CD24","CD27","EGR1","EGF3","CD38","SDC1","MKI67","BCL6","AICDA"), order = T ,cols = colors.FeaturePlot, raster = F)

### relative proportion
Number.HPV_Cluster <- xtabs(~HPV+Cluster,data=TILBs.BCR10X5.combined@meta.data)
Percent.HPV_Cluster <- data.frame(c(Number.HPV_Cluster["HPVN",],Number.HPV_Cluster["HPVP",]),# rownames
                                  as.character(rep(colnames(Number.HPV_Cluster),dim(Number.HPV_Cluster)[1])), #dim(Number.HPV_Cluster)[1]
                                  c(rep("HPVN",dim(Number.HPV_Cluster)[2]),rep("HPVP",dim(Number.HPV_Cluster)[2]))) #dim(Number.HPV_Cluster)[2]

colnames(Percent.HPV_Cluster) <- c("Number","Cluster","HPV")
Percent.HPV_Cluster$Percent <- "Percent"
id <- sort(colnames(Number.HPV_Cluster))
Percent.HPV_Cluster$Cluster = factor(Percent.HPV_Cluster$Cluster, levels =id) # id, C1-C4

ggplot(data = Percent.HPV_Cluster, mapping =aes(x = HPV, y = Number,fill= Cluster)) +
  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
  scale_fill_manual(values = colors.umap)+
  ylab("Relative proportion")+
  #coord_flip()+ # 反转
  theme_classic(12)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black", size=12),
        axis.text.x = element_text(angle =0,size=12,colour="black",hjust = 0.5,vjust=1),
        axis.text.y = element_text(size=12,colour="black"),
        legend.title=element_text(size=12))# 图例名称 字体大小

### HPVP vs. HPVN DEGs analysis
TILBs.HPV <- subset(FindMarkers(TILBs.BCR10X5.combined,ident.1 = "HPVP",group.by = "HPV",logfc.threshold = 0.5),p_val_adj<0.01)
TILBs.HPV <- TILBs.HPV[order(TILBs.HPV$avg_log2FC,decreasing = TRUE),]
TILBs.HPV$id <- rownames(TILBs.HPV)
#write.csv(C1.HPV,"TILBs.HPV.csv")
TILBs.HPV.hi <- subset(TILBs.HPV,avg_log2FC > 0.5)
TILBs.HPV.lo <- subset(TILBs.HPV,avg_log2FC < -0.5)

TILBs.HPV.allgenes <- subset(FindMarkers(TILBs.BCR10X5.combined,ident.1 = "HPVP",group.by = "HPV",logfc.threshold = 0),p_val_adj<2)
TILBs.HPV.allgenes <- TILBs.HPV.allgenes[order(TILBs.HPV.allgenes$avg_log2FC,decreasing = TRUE),]
TILBs.HPV.allgenes$id <- rownames(TILBs.HPV.allgenes)
TILBs.HPV.allgenes$rank <- seq(1,dim(TILBs.HPV.allgenes)[1])
subset(TILBs.HPV.allgenes,avg_log2FC > 0.5 & p_val_adj < 0.01)
subset(TILBs.HPV.allgenes,avg_log2FC < -0.5 & p_val_adj < 0.01)

Ig.Makrers <- c("IGLV3-19","IGLV3-21","IGLV3-1","IGKV3-20","IGHV3-30","MZB1","IGKV1-5","IGKV3-15","CXCR4","SDC1")
#Low.markers <- c("CD37","CXCR4","GAPDH")
ggplot(TILBs.HPV.allgenes,aes(rank,avg_log2FC,label=id))+geom_point(colour='Gray',size=4,shape=1)+
  geom_point(data =  subset(TILBs.HPV.allgenes, avg_log2FC > 0.5 & p_val_adj <0.01),color="red",size=4,shape=1)+
  geom_point(data = TILBs.HPV.allgenes[Ig.Makrers,],color="red",size=4,shape=2)+
  geom_text_repel(data = TILBs.HPV.allgenes[Ig.Makrers,],color="black",size=4)+
  #geom_point(data =  subset(TILBs.HPV.allgenes, avg_log2FC < -0.5 & p_val_adj <0.01),color="blue",size=4,shape=1)+
  #geom_point(data = TILBs.HPV.allgenes[Low.markers,],color="red",size=4,shape=1)+
  #geom_text_repel(data = TILBs.HPV.allgenes[Low.markers,],color="black",size=4)+
  theme_classic(base_size = 16)

### DEGs enrichment
keytypes(org.Hs.eg.db)

Hseg.TILBs.HPV.hi <- bitr(TILBs.HPV.hi$id, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
Hseg.TILBs.HPV.hi.GO <- enrichGO(Hseg.TILBs.HPV.hi$ENTREZID, OrgDb = org.Hs.eg.db, ont='All',pAdjustMethod = 'BH',pvalueCutoff = 0.05, qvalueCutoff = 0.2,keyType = 'ENTREZID',readable = TRUE)
Hseg.TILBs.HPV.hi.GO@result
#write.csv(Hseg.TILBs.HPV.hi.GO@result,"Hseg.TILBs.HPV.hi.GO.csv")
barplot(Hseg.TILBs.HPV.hi.GO, split="ONTOLOGY")+ facet_grid(ONTOLOGY~., scale="free")
Hseg.TILBs.HPV.hi.GO@result$rank <- seq(1,dim(Hseg.TILBs.HPV.hi.GO@result)[1])

GeneRatio <- c()
for (i in Hseg.TILBs.HPV.hi.GO@result$GeneRatio){
  temp1 <- as.numeric(strsplit(i,split = "/")[[1]][1])
  temp2 <- as.numeric(strsplit(i,split = "/")[[1]][2])
  temp <- temp1/temp2
  GeneRatio <- c(GeneRatio,temp)
}

GeneRatio <- round(GeneRatio,2)
Hseg.TILBs.HPV.hi.GO@result$GeneRatio <- GeneRatio
Hseg.TILBs.HPV.hi.GO@result$Description = factor(Hseg.TILBs.HPV.hi.GO@result$Description , levels =  Hseg.TILBs.HPV.hi.GO@result$Description)

meta <- head(Hseg.TILBs.HPV.hi.GO@result,50)
meta <- meta[-8,]
dim(meta)
selection <- c(1,3,4,5,17,23)
ggplot(meta[selection,] ,aes(x=-log10(pvalue),y=Description))+
  labs(title="GO: Biological Process")+
  expand_limits(x=c(8,12))+
  theme_classic(base_size = 20)+
  geom_point(aes(size = Count,colour=GeneRatio), shape = 16)+  scale_size(range = c(5, 10))+scale_color_gradient(low = "#cf0000",high = "#f5c6c6")#color = "Blue

### rds save
#saveRDS(TILBs.BCR10X5.combined,"rds/TILBs.BCR10X5.combinedHPV.PC6.res04.rds")
#saveRDS(Hseg.TILBs.HPV.hi.GO,"rds/Hseg.TILBs.HPV.hi.GO.rds")



### Integrate BCR clone info into meta@data
TILBs.BCR10X5.combined <- combineExpression(BCR.2,TILBs.BCR10X5.combined,
                                 cloneCall="gene+nt", group.by = "sample", proportion = FALSE,
                                 cloneTypes = c(None = 0, Single = 1, Small = 5, Medium = 10,Large = 50, Hyperexpanded = 100),
                                 filterNA = FALSE)

TILBs.BCR10X5.combined@meta.data$BCR_number <- "BCR_number"
TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,cloneType=="Single (0 < X <= 1)")),]$BCR_number <- "BCR=1"
TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,cloneType!="Single (0 < X <= 1)")),]$BCR_number <- "BCR>1"
#TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,Frequency==1)),]$BCR_number <- "BCR=1"
#TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,Frequency>1 & Frequency<10 )),]$BCR_number <- "1<BCR<10"
#TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,Frequency>=10 )),]$BCR_number <- "BCR>=10"
#subset(TILBs.BCR10X5.combined@meta.data, BCR_number=="BCR=1")

CTstrict <- data.frame(xtabs(~CTstrict,TILBs.BCR10X5.combined@meta.data))
CTstrict <- CTstrict[order(CTstrict$Freq,decreasing = TRUE),]

TILBs.BCR10X5.combined@meta.data$TopBCR <- "TopBCR"
TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,CTstrict=="NA_NA_IGLC:LD.7_IGKV3-20")),]$TopBCR <- "BCR-1"
TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,CTstrict=="NA_NA_IGLC:LD.7_IGKV3-15")),]$TopBCR <- "BCR-2"
TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,CTstrict=="NA_NA_IGLC:LD.7_IGKV1D-39")),]$TopBCR <- "BCR-3"
TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,CTstrict=="NA_NA_IGLC:LD.42_IGLV3-21")),]$TopBCR <- "BCR-4"
TILBs.BCR10X5.combined@meta.data[rownames(subset(TILBs.BCR10X5.combined@meta.data,CTstrict=="NA_NA_IGLC:LD.7_IGKV4-1")),]$TopBCR <- "BCR-5"

xtabs(~TopBCR+Cluster,TILBs.BCR10X5.combined@meta.data)
xtabs(~TopBCR,TILBs.BCR10X5.combined@meta.data)

### plot
Idents(TILBs.BCR10X5.combined) <- "Cluster"
DimPlot(TILBs.BCR10X5.combined, reduction = "umap",pt.size = 0.5,order = TRUE,label.size = 4,group.by = "BCR_number",cols = colors.BCR)

Idents(TILBs.BCR10X5.combined) <- "TopBCR"
levels(TILBs.BCR10X5.combined) <- c("BCR-1","BCR-2", "BCR-3","BCR-4","BCR-5","TopBCR")
DimPlot(TILBs.BCR10X5.combined, reduction = "umap",pt.size = 0.2,order = F,label.size = 4,group.by = "TopBCR",cols = colors.umap[c(1,2,3,5,6,4)])

### relative proportion
Number.Cluster_BCRnumber <- xtabs(~TopBCR+Cluster,data=TILBs.BCR10X5.combined@meta.data)[-6,]
Percent.Cluster_BCRnumber <- data.frame(c(Number.Cluster_BCRnumber["BCR-1",],Number.Cluster_BCRnumber["BCR-2",],Number.Cluster_BCRnumber["BCR-3",],Number.Cluster_BCRnumber["BCR-4",],Number.Cluster_BCRnumber["BCR-5",]),# rownames
                                        as.character(rep(colnames(Number.Cluster_BCRnumber),dim(Number.Cluster_BCRnumber)[1])), #dim(Number.Cluster_BCRnumber)[1]
                                        c(rep("BCR-1",dim(Number.Cluster_BCRnumber)[2]),rep("BCR-2",dim(Number.Cluster_BCRnumber)[2]),rep("BCR-3",dim(Number.Cluster_BCRnumber)[2]),
                                          rep("BCR-4",dim(Number.Cluster_BCRnumber)[2]),rep("BCR-5",dim(Number.Cluster_BCRnumber)[2]))) #dim(Number.Cluster_BCRnumber)[2]

colnames(Percent.Cluster_BCRnumber) <- c("Number","Cluster","TopBCR")
#id <- sort(colnames(Number.Cluster_BCRnumber))
#Percent.Cluster_BCRnumber$Cluster = factor(Percent.Cluster_BCRnumber$Cluster, levels =id) #

ggplot(data = Percent.Cluster_BCRnumber, mapping =aes(x = TopBCR, y = Number,fill=Cluster)) +
  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
  scale_fill_manual(values = colors.umap)+
  ylab("Relative proportion")+
  #geom_text(aes(label=number), position=position_stack(vjust=0.5)) +
  #coord_flip()+ # 反转
  #xlab("Patient")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black", size=10),
        axis.text.x = element_text(angle =0,size=10,colour="Black",hjust = 0.5,vjust=0.5),
        axis.text.y = element_text(size=10,colour="Black"),
        legend.title=element_text(size=10))# 图例名称 字体大小

Number.Cluster_BCRnumber <- xtabs(~Cluster+BCR_number,data=TILBs.BCR10X5.combined@meta.data)
Percent.Cluster_BCRnumber <- data.frame(c(Number.Cluster_BCRnumber["C1",],Number.Cluster_BCRnumber["C2",],Number.Cluster_BCRnumber["C3",],Number.Cluster_BCRnumber["C4",],Number.Cluster_BCRnumber["C5",],Number.Cluster_BCRnumber["C6",]),# rownames
                                  as.character(rep(colnames(Number.Cluster_BCRnumber),dim(Number.Cluster_BCRnumber)[1])), #dim(Number.Cluster_BCRnumber)[1]
                                  c(rep("C1",dim(Number.Cluster_BCRnumber)[2]),rep("C2",dim(Number.Cluster_BCRnumber)[2]),rep("C3",dim(Number.Cluster_BCRnumber)[2]),
                                    rep("C4",dim(Number.Cluster_BCRnumber)[2]),rep("C5",dim(Number.Cluster_BCRnumber)[2]),rep("C6",dim(Number.Cluster_BCRnumber)[2]))) #dim(Number.Cluster_BCRnumber)[2]

colnames(Percent.Cluster_BCRnumber) <- c("Number","BCR_number","Cluster")
#id <- sort(colnames(Number.Cluster_BCRnumber))
#Percent.Cluster_BCRnumber$Cluster = factor(Percent.Cluster_BCRnumber$Cluster, levels =id) #

ggplot(data = Percent.Cluster_BCRnumber, mapping =aes(x = Cluster, y = Number,fill= BCR_number)) +
  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
  scale_fill_manual(values = colors.BCR)+
  ylab("Relative proportion")+
  #coord_flip()+ # 反转
  theme_classic(12)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black", size=12),
        axis.text.x = element_text(angle =0,size=12,colour="Black",hjust = 0.5,vjust=1),
        axis.text.y = element_text(size=12,colour="Black"),
        legend.title=element_text(size=10)) # 图例字体大小

Number.Patient_BCRnumber <- xtabs(~Patient+BCR_number,data=TILBs.BCR10X5.combined@meta.data)
Percent.Patient_BCRnumber <- data.frame(c(Number.Patient_BCRnumber["P7",],Number.Patient_BCRnumber["P8",],Number.Patient_BCRnumber["P9",],Number.Patient_BCRnumber["P10",]),# rownames
                                        as.character(rep(colnames(Number.Patient_BCRnumber),dim(Number.Patient_BCRnumber)[1])), #dim(Number.Patient_BCRnumber)[1]
                                        c(rep("P7",dim(Number.Patient_BCRnumber)[2]),rep("P8",dim(Number.Patient_BCRnumber)[2]),rep("P9",dim(Number.Patient_BCRnumber)[2]),
                                          rep("P10",dim(Number.Patient_BCRnumber)[2]))) #dim(Number.Patient_BCRnumber)[2]

colnames(Percent.Patient_BCRnumber) <- c("Number","BCR_number","Patient")
Percent.Patient_BCRnumber$Patient = factor(Percent.Patient_BCRnumber$Patient, levels =c("P7","P8","P9","P10")) #

ggplot(data = Percent.Patient_BCRnumber, mapping =aes(x = Patient, y = Number,fill= BCR_number)) +
  geom_bar(stat= 'identity',position = 'fill',colour= 'black')+
  scale_fill_manual(values = colors.BCR)+
  ylab("Relative proportion")+
  #coord_flip()+ # 反转
  theme_classic(12)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(colour="black", size=12),
        axis.text.x = element_text(angle =0,size=12,colour="Black",hjust = 0.5,vjust=1),
        axis.text.y = element_text(size=12,colour="Black"),
        legend.title=element_text(size=12))# 图例名称 字体大小

### BCR diversity analysis
BCR.meta <- expression2List(TILBs.BCR10X5.combined, split.by = "Cluster")

BCR.meta_P7 <-  expression2List(subset(TILBs.BCR10X5.combined,Patient=="P7"),split.by = "orig.ident")
BCR.meta_P8 <-  expression2List(subset(TILBs.BCR10X5.combined,Patient=="P8"),split.by = "orig.ident")
BCR.meta_P9 <-  expression2List(subset(TILBs.BCR10X5.combined,Patient=="P9"),split.by = "orig.ident")
BCR.meta_P10 <-  expression2List(subset(TILBs.BCR10X5.combined,Patient=="P10"),split.by = "orig.ident")

clonalDiversity(BCR.meta_P7, cloneCall = "gene+nt",group.by = "orig.ident",exportTable = T)
clonalDiversity(BCR.meta_P8, cloneCall = "gene+nt",group.by = "orig.ident",exportTable = T)
clonalDiversity(BCR.meta_P9, cloneCall = "gene+nt",group.by = "orig.ident",exportTable = T)
clonalDiversity(BCR.meta_P10, cloneCall = "gene+nt",group.by = "orig.ident",exportTable = T)

BCR.meta_P7 <- expression2List(subset(TILBs.BCR10X5.combined,Patient=="P7"),split.by = "Cluster")
BCR.meta_P8 <- expression2List(subset(TILBs.BCR10X5.combined,Patient=="P8"),split.by = "Cluster")
BCR.meta_P9 <- expression2List(subset(TILBs.BCR10X5.combined,Patient=="P9"),split.by = "Cluster")
BCR.meta_P10 <- expression2List(subset(TILBs.BCR10X5.combined,Patient=="P10"),split.by = "Cluster")

clonalDiversity(BCR.meta_P7, cloneCall = "gene+nt",group.by = "Cluster",exportTable = T)
clonalDiversity(BCR.meta_P8, cloneCall = "gene+nt",group.by = "Cluster",exportTable = T)
clonalDiversity(BCR.meta_P9, cloneCall = "gene+nt",group.by = "Cluster",exportTable = T)
clonalDiversity(BCR.meta_P10, cloneCall = "gene+nt",group.by = "Cluster",exportTable = T)

clonalDiversity(BCR.meta_P7, cloneCall = "nt",group.by = "Cluster",exportTable = T)
clonalDiversity(BCR.meta_P8, cloneCall = "nt",group.by = "Cluster",exportTable = T)
clonalDiversity(BCR.meta_P9, cloneCall = "nt",group.by = "Cluster",exportTable = T)
clonalDiversity(BCR.meta_P10, cloneCall = "nt",group.by = "Cluster",exportTable = T)

### rds save
#saveRDS(TILBs.BCR10X5.combined,"rds/TILBs.BCR10X5.combinedHPV.PC6.res04.BCRclone.rds")
