library(SeuratObject) # v4.1.3
library(Seurat) # v4.3.0
library(MySeuratWrappers)
library(ggplot2)
library(paletteer)
library(ggrepel)
library(infercnv)
library(phylogram)

### Rscirpt for Figure 3A, 3B, 3G, S6A, S6B, S6C

### define colors panel
colors.FeaturePlot <- c("lightgrey","red")
colors.umap <- paletteer_d("ggsci::default_igv")[c(1,2,3,5)]
colors.Patient <- paletteer_d("ggthemes::Tableau_10")
colors.Tissue <- c("#009966FF","#D60047FF")
#colors.HPV <-  c("#E6BB99","#C37A78")
colors.HPV <-  c("#660099FF","#FF1463FF")

### extract epithelial cells and infer maligant cells from NATs and tumor 10X3 scRNAseq data
PSCC10X3 <- readRDS("rds/PSCC10X3.rmMtRb.combinedPatient.PC10.res01.rds")
Epi10X3 <- subset(PSCC10X3, CellTypes=="Epithelial cells")

Epi10X3.data <- as.matrix(GetAssayData(Epi10X3,slot = "counts"))
Epi10X3.data <- Epi10X3.data[which(rowSums(Epi10X3.data) > 0),] # remove zero genes

Epi10X3.geneOrdering <- read.table("cellranger-6.0.1/refdata-gex-GRCh38-2020-A/chromosome_inferCNV.txt",row.names = 1) #36549 genes
Epi10X3.geneOrdering <- Epi10X3.geneOrdering[rownames(Epi10X3.data),]
Epi10X3.geneOrdering <- na.omit(Epi10X3.geneOrdering)

chr <- c()
for (i in Epi10X3.geneOrdering$V2 ){
  #chr <- c(chr,as.numeric(strsplit(i,"chr")[[1]][2]))
  chr <- c(chr,strsplit(i,"chr")[[1]][2])
}

Epi10X3.geneOrdering$rank <- chr
Epi10X3.geneOrdering <- Epi10X3.geneOrdering[order(Epi10X3.geneOrdering$rank,Epi10X3.geneOrdering$V3),]
Epi10X3.data <- Epi10X3.data[rownames(Epi10X3.geneOrdering),]

Epi10X3.Anno <- data.frame(rownames(Epi10X3@meta.data),Epi10X3@meta.data$Tissue)
colnames(Epi10X3.Anno) <- c("Barcode","Tissue")

#write.table(Epi10X3.Anno,"Epi10X3.Anno.txt",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE) #  Manual sort by chr
#write.table(Epi10X3.geneOrdering,"Epi10X3.geneOrdering.txt",sep = "\t",quote = FALSE,row.names = TRUE,col.names = FALSE)

### inferCNV analysis
chrOrdering <- read.table("Epi10X3.geneOrdering.txt",row.names = 1)
Epi10X3.data <- Epi10X3.data[rownames(chrOrdering),]

Epit10X3.infercnv_obj = CreateInfercnvObject(raw_counts_matrix=Epi10X3.data,
                                         annotations_file="Epi10X3.Anno.txt",
                                         delim="\t",
                                         gene_order_file="Epi10X3.geneOrdering.txt",
                                         ref_group_names=c("NATs"))

Epit10X3.infercnv_obj = infercnv::run(Epit10X3.infercnv_obj,
                                  cutoff=0.1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                                  out_dir="inferCNV_NATs_vs_Tumor",  # Auto create document
                                  cluster_by_groups=T,   # cluster
                                  denoise=T, #去噪
                                  HMM=T) # 是否基于HMM预测CNV



### extract maligant epithelial cells and reanalysis from tumor 10X3 scRNAseq data
PSCC10X3 <- readRDS("rds/PSCC10X3.rmMtRb.combinedPatient.PC10.res01.rds")
Epi10X3 <- subset(PSCC10X3, CellTypes=="Epithelial cells")
EpiTumor10X3<- subset(Epi10X3, Tissue=="Tumor")
#xtabs(~Patient, EpiTumor10X3@meta.data)

Epi.data <- GetAssayData(object = EpiTumor10X3, slot = "counts")
EpiTumor10X3 <- CreateSeuratObject(counts = Epi.data, min.cells = 10, min.features = 200, project = "EpiTumor10X3")
EpiTumor10X3@meta.data <- PSCC10X3@meta.data[rownames(EpiTumor10X3@meta.data),]

###  Integrate data
Idents(EpiTumor10X3) <- "Patient"
new.cluster.ids <- c("T1","T2","T3","T2","T2","T6") # Integrate P2,P4,P5 into T2 to find anchors due to less cells
names(new.cluster.ids) <- levels(EpiTumor10X3)
EpiTumor10X3 <- RenameIdents(EpiTumor10X3,new.cluster.ids)
EpiTumor10X3@meta.data["Trank"] <- data.frame(Idents(EpiTumor10X3))[1]
EpiTumor10X3@meta.data <- EpiTumor10X3@meta.data[c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Patient","Tissue","HPV","Groups","Library","Trank")]

batch.list <- SplitObject(EpiTumor10X3,split.by = "Trank")
batch.list <- lapply(X = batch.list, FUN = function(x) {
  x <- NormalizeData(x,normalization.method = "LogNormalize", scale.factor = 10000)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures =2000)
})
features <- SelectIntegrationFeatures(object.list = batch.list)
anchors <- FindIntegrationAnchors(object.list = batch.list, anchor.features = features)
EpiTumor10X3.combined <- IntegrateData(anchorset = anchors)

### UMAP analysis
DefaultAssay(EpiTumor10X3.combined) <- "integrated"
EpiTumor10X3.combined <- ScaleData(EpiTumor10X3.combined,features = rownames(EpiTumor10X3.combined))
EpiTumor10X3.combined <- RunPCA(EpiTumor10X3.combined,features = VariableFeatures(EpiTumor10X3.combined))

EpiTumor10X3.combined <- JackStraw(EpiTumor10X3.combined,dims = 50)
EpiTumor10X3.combined <- ScoreJackStraw(EpiTumor10X3.combined , dims = 1:50)
JackStrawPlot(EpiTumor10X3.combined, dims = 1:20)
ElbowPlot(EpiTumor10X3.combined)

EpiTumor10X3@meta.data <- EpiTumor10X3@meta.data[c("orig.ident","nCount_RNA","nFeature_RNA","percent.mt","Patient","Tissue","HPV","Groups","Library")]
#EpiTumor10X3.combined <- RunTSNE(EpiTumor10X3.combined, dims = 1:8)
EpiTumor10X3.combined <- RunUMAP(EpiTumor10X3.combined, dims = 1:8)
EpiTumor10X3.combined <- FindNeighbors(EpiTumor10X3.combined,reduction = "pca",dims = 1:8)
EpiTumor10X3.combined <- FindClusters(EpiTumor10X3.combined, resolution =0.1)

### cell labels
Idents(EpiTumor10X3.combined) <- "integrated_snn_res.0.1"
new.cluster.ids <- c("C2","C1","C3","C4")
names(new.cluster.ids) <- levels(EpiTumor10X3.combined)
EpiTumor10X3.combined <- RenameIdents(EpiTumor10X3.combined,new.cluster.ids)
EpiTumor10X3.combined@meta.data["Cluster"] <- data.frame(Idents(EpiTumor10X3.combined))[1]
levels(EpiTumor10X3.combined) <- c("C1","C2","C3","C4")

### Find Cluster Markers
Idents(Immu10X3.combined) <- "Cluster"
DefaultAssay(EpiTumor10X3.combined) <- "RNA"
subtype.markers <- FindAllMarkers(EpiTumor10X3.combined,logfc.threshold = 0.5,only.pos = TRUE) # Rank Cluster#
subset(subtype.markers,cluster=="C2")
subset(subtype.markers,gene=="TOP2A")

### plot
DimPlot(EpiTumor10X3.combined,reduction = "umap",order=T, raster = T, label = F, cols = colors.umap,repel = T,label.size = 4,)
DimPlot(EpiTumor10X3.combined,reduction = "umap",order=T, raster = T,label = F,cols = colors.HPV, repel = T,group.by = "HPV",label.size = 4)
DimPlot(EpiTumor10X3.combined,reduction = "umap",order=T, raster = T,label = F,cols = colors.Patient, repel = T,group.by = "Patient",label.size = 4)

DefaultAssay(EpiTumor10X3.combined) <- "RNA"
subtype.markers.top <-  subtype.markers %>% group_by(cluster) %>% top_n(20, avg_log2FC)
levels(EpiTumor10X3.combined) <- c("C1","C2","C3","C4")
DotPlot(EpiTumor10X3.combined,features = subtype.markers.top$gene,cols = colors.FeaturePlot )+RotatedAxis()# + coord_flip()

MySeuratWrappers::VlnPlot(EpiTumor10X3.combined,features = c("PIGR","LCN2","PSCA","SELE","ACKR1","TRIM31","TMC5"),pt.size =0,stack=T,x.lab = "",y.lab = "",cols =colors.umap)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

### relative proportion
Number.HPV_Cluster <- xtabs(~HPV+Cluster,data=EpiTumor10X3.combined@meta.data)
Percent.HPV_Cluster <- data.frame(c(Number.HPV_Cluster["HPVN",],Number.HPV_Cluster["HPVP",]),# rownames
                                  as.character(rep(colnames(Number.HPV_Cluster),dim(Number.HPV_Cluster)[1])), #dim(Number.HPV_Cluster)[1]
                                  c(rep("HPVN",dim(Number.HPV_Cluster)[2]),rep("HPVP",dim(Number.HPV_Cluster)[2]))) #dim(Number.HPV_Cluster)[2]

colnames(Percent.HPV_Cluster) <- c("Number","Cluster","HPV")
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
        axis.text.x = element_text(angle =45,size=12,colour="Black",hjust = 1,vjust=1),
        axis.text.y = element_text(size=12,colour="Black"),
        legend.title=element_text(size=12))# 图例名称 字体大小

### Cluster DEGs analysis
C3V2.fc0 <- FindMarkers(EpiTumor10X3.combined,ident.1 = "C3",ident.2 = "C2",logfc.threshold = 0)
C3V2.fc0["gene"] <- rownames(C3V2.fc0)
C3V2 <- subset(C3V2.fc0,p_val_adj<0.01 & avg_log2FC>1) # 172
C3V2 <- C3V2[order(C3V2$avg_log2FC,decreasing = TRUE),]
C3V2["PIGR",]
C3V2.up <- (C3V2 %>% top_n(10, avg_log2FC))$gene
C3V2.down <- (C3V2 %>% top_n(-10, avg_log2FC))$gene

Idents(EpiTumor10X3.combined) <- "Cluster"
C4V2.fc0 <- FindMarkers(EpiTumor10X3.combined,ident.1 = "C4",ident.2 = "C2",logfc.threshold = 0)
C4V2.fc0["gene"] <- rownames(C4V2.fc0)
C4V2 <- subset(C4V2.fc0,p_val_adj<0.01 & avg_log2FC>1) # 289
C4V2 <- C4V2[order(C4V2$avg_log2FC,decreasing = TRUE),]
C4V2["PIGR",]
C4V2.up <- (C4V2 %>% top_n(10, avg_log2FC))$gene
C4V2.down <- (C4V2 %>% top_n(-10, avg_log2FC))$gene

intersect(rownames(C4V2),rownames(C3V2)) # 107
setdiff(rownames(C4V2),rownames(C3V2)) # 182
setdiff(rownames(C3V2),rownames(C4V2)) # 65

C34V2.fc0 <- FindMarkers(EpiTumor10X3.combined,ident.1 = c("C3","C4"),ident.2 = "C2",logfc.threshold = 0)
C34V2.fc0["gene"] <- rownames(C34V2.fc0)
#C34V2 <- FindMarkers(Epit.combined,ident.1 = c("2","3"),ident.2 = "1",logfc.threshold = 0.5,only.pos = TRUE)
C34V2 <- subset(C34V2.fc0,p_val_adj<0.01 & abs(avg_log2FC)>1)
C34V2 <- C34V2[order(C34V2$avg_log2FC,decreasing = TRUE),]
C34V2["PIGR",]
C34V2.up <- (C34V2 %>% top_n(30, avg_log2FC))$gene
C34V2.down <- (C34V2 %>% top_n(-30, avg_log2FC))$gene

ggplot(data = C3V2.fc0, mapping = aes(x = 100*(pct.1-pct.2), y = avg_log2FC,label=gene))+ geom_point(color="Gray",size=0.5,)+
  geom_point(data = subset(C3V2.fc0, C3V2.fc0$avg_log2FC> 1),color="red",size=1,)+ #size是圆圈大小
  geom_point(data = subset(C3V2.fc0, C3V2.fc0$avg_log2FC< -1),color="blue",size=1,)+ #size是圆圈大小
  #geom_point(data = C3V2[C3V2.up,],color="red",size=1,)+
  geom_point(data = C3V2["PIGR",],color="#401c44",size=3,)+
  #geom_text_repel(data=C3V2[C3V2.up,],size=3,color="black",)+ # size是字体大小
  geom_text_repel(data=C3V2["PIGR",],size=4,color="black",)+ # size是字体大小
  #geom_text_repel(data=C3V2[C3V2.down,],size=3,color="black",)+ # size是字体大小
  #scale_y_continuous(limits = c(-4,3))+
  xlab("Δ Difference (%)")+
  ylab("Log-Fold Change")+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0,size=12,colour="Black",hjust = 1),
    axis.title.x = element_text(colour="black", size=14), #坐标轴字体大小
    axis.title.y = element_text(colour="black", size=14),
    axis.text.y = element_text(size=14,colour="Black"),
    legend.title=element_text(size=12),# 图例名称 字体大小
    legend.text=element_text(size=12))# 图例备注 字体大小

ggplot(data = C4V2.fc0, mapping = aes(x = 100*(pct.1-pct.2), y = avg_log2FC,label=gene))+ geom_point(color="Gray",size=0.5,)+
  geom_point(data = subset(C4V2.fc0, C4V2.fc0$avg_log2FC> 1),color="red",size=1,)+ #size是圆圈大小
  geom_point(data = subset(C4V2.fc0, C4V2.fc0$avg_log2FC< -1),color="blue",size=1,)+ #size是圆圈大小
  #geom_point(data = C4V2[C4V2.up,],color="red",size=1,)+
  geom_point(data = C4V2["PIGR",],color="#401c44",size=3,)+
  #geom_text_repel(data=C4V2[C4V2.up,],size=3,color="black",)+ # size是字体大小
  geom_text_repel(data=C4V2["PIGR",],size=4,color="black",)+ # size是字体大小
  #geom_text_repel(data=C4V2[C4V2.down,],size=3,color="black",)+ # size是字体大小
  #scale_y_continuous(limits = c(-4,3))+
  xlab("Δ Difference (%)")+
  ylab("Log-Fold Change")+
  theme_classic()+
  theme(#axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 0,size=12,colour="Black",hjust = 1),
    axis.title.x = element_text(colour="black", size=14), #坐标轴字体大小
    axis.title.y = element_text(colour="black", size=14),
    axis.text.y = element_text(size=14,colour="Black"),
    legend.title=element_text(size=12),# 图例名称 字体大小
    legend.text=element_text(size=12))# 图例备注 字体大小

### rds save
#saveRDS(EpiTumor10X3.combined,"rds/EpiTumor10X3.combined.rmMtRb.combinedTrank.PC8.res01.rds")
#saveRDS(subtype.markers,"rds/subtype.markers.Cluster.EpiTumor10X3.combined.rmMtRb.combinedTrank.PC8.res01.rds")
