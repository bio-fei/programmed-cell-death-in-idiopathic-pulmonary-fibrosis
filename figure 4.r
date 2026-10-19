#A
library(dplyr)
library(magrittr)
library(data.table)
library(Seurat)
library(Matrix)
library(plyr)
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
library(ggsci)
library(future)
library(scattermore)
library(clustree)
library(patchwork)
library(pheatmap)
library(harmony)
set.seed(7889)
plan(multisession, workers = 4)
options(future.globals.maxSize = 450*1000 * 1024^2)
load("FGSE128033_scRNA.RData")
data<- data_GSE128033
Idents(data)<-data$orig.ident
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
data[["percent.ribo"]] <-PercentageFeatureSet(data, "^RP[SL]")
data[["percent.hb"]] <- PercentageFeatureSet(data, "^HB[AB]")
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.ribo","percent.hb"), ncol = 5)
selected_c <- WhichCells(data, expression = nFeature_RNA > 200 & nFeature_RNA <=6000 &  percent.mt <20) # 
gene5    <- rownames(data)[Matrix::rowSums(data) > 4]
hb.genes <-as.data.frame(table(substr(rownames(data)[grep(c("^RP"),rownames(data))],1,4)))
hb.genes <- c(rownames(data)[grep("^RP[SL]",rownames(data))],
              rownames(data)[grep(c("^RP11"),rownames(data))],
              rownames(data)[grep(c("^RP1-"),rownames(data))],
              rownames(data)[grep(c("^RP3-"),rownames(data))],
              rownames(data)[grep(c("^RP4-"),rownames(data))],
              rownames(data)[grep(c("^RP5-"),rownames(data))],
              rownames(data)[grep(c("^RP6-"),rownames(data))],
              rownames(data)[grep(c("^RP13"),rownames(data))])
selected_f <- gene5
data<- subset(data, cells = selected_c, features = selected_f)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data, vars.to.regress = c("percent.mt", "percent.ribo","percent.hb"))
data <- RunPCA(data, features = VariableFeatures(data))
ElbowPlot(data,ndims = 50)
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
cumu <- cumsum(pct)  
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
data <- FindNeighbors(data, dims = 1:pc.use, reduction = "pca")
data <- RunUMAP(data, dims = 1:pc.use, reduction = "pca", reduction.name = "umap.unintegrated")
data1 <- RunHarmony(data,reduction = "pca",group.by.vars = "patient",reduction.save = "harmony")
pct <- data1[["harmony"]]@stdev / sum(data1[["harmony"]]@stdev) * 100
cumu <- cumsum(pct)  
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)

data1 <- FindNeighbors(data1, reduction = "harmony", dims = 1:pc.use)
data1 <- RunUMAP(data1, reduction = "harmony", dims = 1:pc.use, reduction.name = "umap.harmony")
data1 <- FindClusters( object = data1,resolution = c(seq(1,3,0.1)))
i<-0.1
data1 <- FindClusters(data1, resolution =i)
Idents(data1)<-data1$seurat_clusters
marker_gene<-c("EPCAM", # epithelial
               "CLDN5", # endothelial
               "COL1A2",# mesenchyme
               "PTPRC", # immune CD45—official gene name
               "FCER1G") #myeloid
mycolor <- c('white', '#ac5e74')
DotPlot(object = data1,features = marker_gene,cols= mycolor) + RotatedAxis()+coord_flip()

data1@meta.data$cell_type<-""
data1@meta.data$cell_type[which(data1@meta.data$seurat_clusters %in% c("3","10"))]<-"Endothelial"
data1@meta.data$cell_type[which(data1@meta.data$seurat_clusters %in% c("4","6"))]<-"Epithelial"
data1@meta.data$cell_type[which(data1@meta.data$seurat_clusters %in% c("1","8"))]<-"Lymphoid"
data1@meta.data$cell_type[which(data1@meta.data$seurat_clusters %in% c("0","2","7","9"))]<-"Myeloid"
data1@meta.data$cell_type[which(data1@meta.data$seurat_clusters %in% c("5"))]<-"Mesenchymal"
pdf(height = 8, width = 10, file = "9-umap-cell-harmony.pdf")
DimPlot(data1, reduction = "umap.harmony", group.by ="cell_type", label = TRUE, 
        pt.size = 0.2, label.size = 4, label.box = F) +scale_color_npg()
dev.off()
DimPlot(data1, reduction = "umap.unintegrated",  group.by ="cell_type", label = TRUE, 
        pt.size = 0.2, label.size = 4, label.box = F) +scale_color_npg()



#B
Idents(data1)<-data1$cell_type
markers <- FindAllMarkers(data1,pos=T)
library(scRNAtoolVis)
pal<-c(pal_d3("category20")(20)[1:14],pal_locuszoom()(7)[6],pal_d3("category20")(20)[15:20])
p<-jjVolcano(diffData = markers,
             log2FC.cutoff = 1, 
             size  = 3.5, 
             fontface = 'italic', 
             #aesCol = c('purple','orange'), 
             tile.col = pal, 
             #col.type = "adjustP", 
             topGeneN = 5 )


#C
i<-2.1
data1 <- FindClusters(data1, resolution =i)
DimPlot(data1, reduction = "umap.harmony", label = TRUE, pt.size =0.2, label.size = 5, label.box = F) 

marker_gene<-c("GNLY","NKG7","CD3D","THEMIS","CD79A","MS4A1",#Lymphoid
               "CD1E","CSF2RA","S100A12","FCN1", "CD68","MARCO","TPSAB1","TPSB2",# Myeloid
               "MUC5B","ANKRD36C", "RSPH1","FAM183A","AGER","RTKN2","NAPSA","SFTPC","SCGB3A2", "KRT5","KRT17", #Epithelial
               "KCNK3","PDGFRB","ACTA2","TAGLN", "COL1A1","PDGFRA",  #Mesenchymal
               "FCN3","RGCC","CCL21","MMRN1","CLDN5","DKK2","VWF","ACKR1") #Endothelial

marker_gene<-c("MUC5B","ANKRD36C",  #goblet	"RIMS1",
               "KRT5","KRT17",      #Basal "TP63",
               "NAPSA","SFTPC",     #AT2	"HHIP","SFTPD",
               "AGER","RTKN2",      #AT1	 "CLIC5",
               "SCGB3A2",           #club	
               "RSPH1","FAM183A",   #Ciliated	 "FOXJ1","CCDC78",
               "RGS5","PDGFRB",   #Pericyte	
               "ACTA2","TAGLN",    #Smooth Muscle	
               "COL1A1","PDGFRA",  #fibroblast	
               "WT1","CFB",       #mesothelial	
               "FCN3","BTNL9",        #Capillary	
               "SEMA3G","DKK2",         #Artery	
               "VWF","ACKR1",        #Vein	
               "CCL21","MMRN1",        #lymphatic	
               "GNLY","NKG7",                          #NK	
               "CD3D","THEMIS",                         #T 
               "CD79A","MS4A1",                         #B/Plasma	
               "CD1E","CSF2RA",          #Dendritic	
               "S100A12","FCN1",        #Monocyte	
               "CD68","MARCO",          #Macrophage 	
               "TPSAB1","TPSB2")                  #Mast	
mycolor <- c('#F0F8FF','#ac5e74')
DotPlot(object = data1,features = marker_gene,cols= mycolor) + RotatedAxis()+coord_flip()

data1@meta.data$cell_type1<-"nocell"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("12"))]<-"Goblet"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("22"))]<-"Basal"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("13"))]<-"AT2"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("29"))]<-"AT1"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("24"))]<-"Club"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("15","36"))]<-"Ciliated"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("35"))]<-"Pericyte"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("34"))]<-"Smooth_Muscle"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("8","20","28"))]<-"Fibroblast"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("40"))]<-"Mesothelial"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("23","18"))]<-"Capillary"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("32"))]<-"Artery"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("4","39"))]<-"Vein"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("31"))]<-"Lymphatic"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("14","21"))]<-"NK"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("9","3","25","26"))]<-"T"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("30","37"))]<-"B"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("19"))]<-"Dendritic"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("6"))]<-"Monocyte"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("0","1","2","5","7","10","16","17","27","33","38"))]<-"Macrophage"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("11"))]<-"Mast"
pdf(height = 8, width = 10, file = "11-umap-cell20-harmony.pdf")
DimPlot(data1, reduction = "umap.harmony", group.by ="cell_type1", label = TRUE, 
        pt.size = 0.2, label.size = 4, label.box = F,
        cols=c(pal_d3("category20")(20)[1:14],pal_locuszoom()(7)[6],pal_d3("category20")(20)[15:20]))
dev.off()
DimPlot(data1, reduction = "umap.unintegrated", group.by ="cell_type1", label = TRUE, 
        pt.size = 0.2, label.size = 4, label.box = F,
        cols=c(pal_d3("category20")(20)[1:14],pal_locuszoom()(7)[6],pal_d3("category20")(20)[15:20]))



#D
Idents(data1)<-data1$cell_type1
new_order <- c("NK", "T", "B","Dendritic", "Monocyte", "Macrophage","Mast",
               "Goblet", "Ciliated", "AT1","AT2", "Club", "Basal",
               "Pericyte", "Smooth_Muscle","Fibroblast","Mesothelial","Capillary", "Lymphatic", "Artery", "Vein")
data1@active.ident<- factor(data1@active.ident, levels = new_order)
list_genes=list(Lymphoid= c("GNLY","NKG7","CD3D","THEMIS","CD79A","MS4A1"),
                                  Myeloid= c("CD1E","CSF2RA","S100A12","FCN1", "CD68","MARCO","TPSAB1","TPSB2"),
                                  Epithelial=c("MUC5B","ANKRD36C", "RSPH1","FOXJ1","AGER","RTKN2","NAPSA","SFTPC",
                                               "SCGB3A2", "KRT5","KRT17"), 
                                  Mesenchymal= c("RGS5","PDGFRB","ACTA2","TAGLN", "COL1A1","PDGFRA","WT1","CFB"), 
                                  Endothelial= c("FCN3","BTNL9","CCL21","MMRN1","SEMA3G","DKK2","VWF","ACKR1"))
pdf(height =6, width = 12, file = "12-celltype-markers_dot.pdf")
DotPlot(object = data1,features =list_genes,col.min = 0) +xlab("")+ylab("cell type")+
  RotatedAxis()+scale_color_gradientn(colours = colorRampPalette(colors = c("white","#FFCC33"))(100))+ #颜色
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle = 90)) #轴标签
dev.off()



#E
Metadata<-data1@meta.data
Metadata$status_1[which(Metadata$status_1 %in% c("Control_LOW","Control_UP"))]<-"Control"
Metadata$status_1[which(Metadata$status_1 %in% c("IPF_LOW"))]<-"IPF_LOW"
Metadata$status_1[which(Metadata$status_1 %in% c("IPF_UP"))]<-"IPF_UP"
Metadata$number<-1
Metadata$patient<-factor(Metadata$patient,levels = c(unique(Metadata$patient[which(Metadata$status_1 %in% "Control")]),
                                                     unique(Metadata$patient[which(Metadata$status_1 %in% "IPF_UP")]),
                                                     unique(Metadata$patient[which(Metadata$status_1 %in% "IPF_LOW")])))
pdf("14-celltype5_proportion.pdf",width = 8,height = 6)
ggplot(data=Metadata,aes(x=patient,fill=factor(cell_type)))+geom_bar(position='fill')+
  scale_fill_manual(values = pal_npg()(5)) +
  theme( axis.text.x = element_text(size = 12,angle = 45,hjust = 1)) +
  theme(panel.background=element_rect(fill="white",colour="white"), # 填充框内白色背景和白色的网格线 
        axis.line=element_line(colour="black"))# x,y轴颜色，粗细
dev.off()
pdf("14-celltype20_proportion.pdf",width = 8,height = 6)
ggplot(data=Metadata,aes(x=patient,fill=factor(cell_type1)))+geom_bar(position='fill')+
  scale_fill_manual(values =c(pal_d3("category20")(20)[1:14],pal_locuszoom()(7)[6],pal_d3("category20")(20)[15:20])) +
  theme( axis.text.x = element_text(size = 12,angle = 45,hjust = 1)) +
  theme(panel.background=element_rect(fill="white",colour="white"), # 填充框内白色背景和白色的网格线 
        axis.line=element_line(colour="black"))# x,y轴颜色，粗细
dev.off()


#F
library("sscVis")
source("do.tissueDist.R")
Metadata$meta.cluster<-Metadata$cell_type
Metadata$loc<-Metadata$status_1
OR.CD8.list <- do.tissueDist(cellInfo.tb=Metadata,
                             out.prefix="15-celltype5-OR-2",
                             pdf.width=4,pdf.height=6,verbose=1)
Metadata$meta.cluster<-Metadata$cell_type
Metadata$loc<-Metadata$status
OR.CD8.list <- do.tissueDist(cellInfo.tb=Metadata,
                             out.prefix="15-celltype5-OR-1",
                             pdf.width=4,pdf.height=6,verbose=1)
Metadata$meta.cluster<-Metadata$cell_type1
Metadata$loc<-Metadata$status_1
OR.CD8.list <- do.tissueDist(cellInfo.tb=Metadata,
                             out.prefix="15-celltype20-OR-2",
                             pdf.width=4,pdf.height=6,verbose=1)
Metadata$meta.cluster<-Metadata$cell_type1
Metadata$loc<-Metadata$status
OR.CD8.list <- do.tissueDist(cellInfo.tb=Metadata,
                             out.prefix="15-celltype20-OR-1",
                             pdf.width=4,pdf.height=6,verbose=1)


saveRDS(data1, file = "clustered_data_sub_celltype.rds")
