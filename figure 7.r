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
library(SingleR)
set.seed(7889)
plan(multisession, workers = 4)
options(future.globals.maxSize = 450*1000 * 1024^2)
cell_type_cols<-c(pal_simpsons()(16)[1:10],pal_igv("default")(41)[2:30])
load("GSE264278_scRNA.RData")
data<- data_GSE264278
Idents(data)<-data$orig.ident
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^mt-")
data[["percent.ribo"]] <-PercentageFeatureSet(data, "^Rp[sl]")
data[["percent.hb"]] <- PercentageFeatureSet(data, "^Hb[ab]")
selected_c <- WhichCells(data, expression = percent.mt <20) # 
hb.genes <-as.data.frame(table(substr(rownames(data)[grep(c("^Rp"),rownames(data))],1,4)))
hb.genes <- c(rownames(data)[grep("^Rp[sl]",rownames(data))],
              rownames(data)[grep(c("^Rp11"),rownames(data))],
              rownames(data)[grep(c("^Rp1-"),rownames(data))],
              rownames(data)[grep(c("^Rp3-"),rownames(data))],
              rownames(data)[grep(c("^Rp4-"),rownames(data))],
              rownames(data)[grep(c("^Rp5-"),rownames(data))],
              rownames(data)[grep(c("^Rp6-"),rownames(data))],
              rownames(data)[grep(c("^Rp13"),rownames(data))])
rb.genes <- rownames(data)[grep("^Hb[ab]",rownames(data))]
mt.genes <- rownames(data)[grep("^mt-",rownames(data))]
selected_f <- setdiff(rownames(data),c(hb.genes,rb.genes,mt.genes))
data<- subset(data, cells = selected_c,features =selected_f )
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
data <- ScaleData(data)#vars.to.regress = c("percent.mt", "percent.ribo","percent.hb")
data <- RunPCA(data, features = VariableFeatures(data))
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100
cumu <- cumsum(pct)  #20
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
data <- FindNeighbors(data, dims = 1:pc.use, reduction = "pca")
data1 <- RunHarmony(data,reduction = "pca",group.by.vars = "patient",reduction.save = "harmony")
data1 <- FindNeighbors(data1, reduction = "harmony", dims = 1:pc.use)
data1 <- RunUMAP(data1, reduction = "harmony", dims = 1:pc.use, reduction.name = "UMAP.harmony")
data1 <- FindClusters(data1, resolution =0.5)
mycolor <- c('#F0F8FF','#ac5e74')#设置颜色
marker_gene<-c("Krt8",#"MUC5B","ANKRD36C",  #goblet	"RIMS1",
                  "Krt5",    #Basal "TP63",
                  "Sftpb","Sftpc",     #AT2	"HHIP","SFTPD",
                  "Pdpn","Aqp5",      #AT1	 "CLIC5",
                  "Scgb1a1","Scgb3a2",           #club	
                  "Lrrc23","Foxj1",   #Ciliated	 "FOXJ1","CCDC78",
                  "Pecam1","Rgs5",   #Pericyte	
                  "Acta2","Myl9",    #Smooth Muscle	
                  "Col1a1",  #fibroblast	
                  "Gpihbp1",        #Capillary	
                  "Ptgs1",        #Vein	
                  "Lyve1",        #lymphatic	
                  "Ptprc","Gzma",         #NK	
                  "Cd3e","Cd3g",          #T 
                  "Cd19","Cd79a",          #B/Plasma	
                  "Fcer1g","Cd74","Cd83",          #Dendritic	
                  "Ccr2","Lyz2" ,       #Monocyte	
                  "Cd68",          #Macrophage 	
                  "Fcer1a") #Mast	
DotPlot(object = data1,features = marker_gene,cols= mycolor) + RotatedAxis()+coord_flip()
data1@meta.data$cell_type1<-"nocell"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("9"))]<-"AT2"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("12"))]<-"AT1"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("22"))]<-"Club"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("18"))]<-"Ciliated"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("14","15"))]<-"Mural cells"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("2","13"))]<-"Fibroblast"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("1","17","24"))]<-"vascular endothelial cell"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("19"))]<-"Lymphatic"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("3","7","10"))]<-"NK/T"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("0","23","26"))]<-"B"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("11"))]<-"Dendritic"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("6"))]<-"Monocyte"
data1@meta.data$cell_type1[which(data1@meta.data$seurat_clusters %in% c("25","4","5","16"))]<-"Macrophage"
Idents(data1)<-data1$cell_type1
selected_c <- WhichCells(data1, idents=setdiff(unique(data1$cell_type1),"nocell")) 
data<- subset(data1, cells = selected_c)
pdf(height = 8, width = 10, file = "GSE264278-3-subcell.pdf")
DimPlot(data, reduction = "UMAP.harmony", group.by ="cell_type1", label = TRUE, 
        pt.size = 1, label.size = 4, label.box = F,
        cols=c(pal_simpsons()(16)[1:10],pal_igv("default")(41)[2:30]))
dev.off()

pdf("GSE264278-4-CD24-DotPlot.pdf",width = 10,height = 9)
DotPlot(object = data,features = c("Krt8","Krt18","Cd24a","Lrrc23","Foxj1"),cols= mycolor) + RotatedAxis()+coord_flip()
dev.off()
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])

coef$mmu_gene<-c("Cdh3","Mmp11","Mmp13","Cdh2","Prom2","Il13ra2","Adra1a","Sfrp2",
"Cd24a","Twist1","Cldn1","Pthlh","Fzd5","Eda2r","Cck","Nfatc4","Faim","Cdcp1","Uvrag")
expr_matrix <- GetAssayData(data1, assay = "RNA", slot = "data")  # 使用 slot="data"（log-normalized）或 "counts"
intersect(coef$mmu_gene, rownames(expr_matrix))
expr_matrix<-expr_matrix[intersect(coef$mmu_gene, rownames(expr_matrix)),]
rownames(expr_matrix)[1:5]
rownames(coef)<-coef$mmu_gene
coef1<-coef[intersect(coef$mmu_gene, rownames(expr_matrix)),]
expr_matrix<-as.data.frame(t(as.matrix(expr_matrix)))
for (ii in 1:nrow(expr_matrix)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*expr_matrix[ii,jj]
  }
  expr_matrix[ii,20]<-pcd_score
}
colnames(expr_matrix)[20]<- "PCD_scores"
data1$PCD_scores <- expr_matrix$PCD_scores
p<-ggplot(data1@meta.data,aes(x=cell_type1,y=PCD_scores,fill = cell_type1))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+ 
  theme(legend.position = 'none')+xlab('')+ylab('PCD score')+labs(fill='')+
  scale_color_npg()
ggsave(p,file="GSE264278-7-celltype-PCD_scores.pdf",width = 7,height = 4) 
p<-ggplot(data1@meta.data,aes(x=cell_type1,y=PCD_scores,fill=status))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+#scale_fill_manual(values = c("#66CCFF","#FF6600"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", " ")),label = "p.signif")+
  theme(legend.position = 'top')+xlab('')+ylab('PCD score')+labs(fill='')
ggsave(p,file="GSE264278-8-celltype-status-PCD_scores.pdf",width = 10,height = 4) 


