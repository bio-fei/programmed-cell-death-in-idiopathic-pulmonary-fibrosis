#A
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
data1<-readRDS("clustered_data_sub_celltype.rds")
expr_matrix <- GetAssayData(data1, assay = "RNA", slot = "data")  
intersect(coef$gene, rownames(expr_matrix))
expr_matrix<-expr_matrix[intersect(coef$gene, rownames(expr_matrix)),]
rownames(expr_matrix)[1:5]
coef1<-coef[intersect(coef$gene, rownames(expr_matrix)),]

expr_matrix<-as.data.frame(t(as.matrix(expr_matrix)))
for (ii in 1:nrow(expr_matrix)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*expr_matrix[ii,jj]
  }
  expr_matrix[ii,20]<-pcd_score
}
colnames(expr_matrix)[20]<- "PCD_scores"
sum(expr_matrix$PCD_scores>0)

data1$PCD_scores <- expr_matrix$PCD_scores
pdf(height = 8, width = 8, file = "18-PCD_scores_UMAP.pdf")
FeaturePlot(data1, features = 'PCD_scores',cols = c('#FFFAFA',"#CC0000"), pt.size = 0.2,reduction = "umap.harmony")
dev.off()




#B
my_compare<-list(c('Control','IPF'))
p<-ggplot(Metadata,aes(x=status,y=PCD_scores,fill=status))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF6600"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
 stat_compare_means(comparisons=my_compare,method = "t.test")+
  theme(legend.position = 'top')+xlab('')+ylab('PCD score')+labs(fill='')
ggsave(p,file="19-status-PCD_scores.pdf",width = 7,height = 4)
my_compare<-list(c('Control','IPF_UP'), c('IPF_UP','IPF_LOW'),c('Control','IPF_LOW'))
p<-ggplot(Metadata,aes(x=status_1,y=PCD_scores,fill=status_1))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF9900","#CC3300"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(comparisons=my_compare,method = "t.test")+
  theme(legend.position = 'top')+xlab('')+ylab('PCD score')+labs(fill='')
ggsave(p,file="19-status_1-PCD_scores.pdf",width = 7,height = 4) 




#C
seurat_obj <- Load10X_Spatial(data.dir = "E-MTAB-14121\\V10T31-051-D1")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, "^MT-")
seurat_obj[["percent.ribo"]] <-PercentageFeatureSet(seurat_obj, "^RP[SL]")
seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, "^HB[AB]")
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj,vars.to.regress = "percent.mt")   
seurat_obj <- RunPCA(seurat_obj)
pct <-seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
cumu <- cumsum(pct) 
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:pc.use)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pc.use)
expr_matrix <- GetAssayData(seurat_obj, assay = "Spatial", slot = "data")  # 使用 slot="data"（log-normalized）或 "counts"
intersect(coef$gene, rownames(expr_matrix))
expr_matrix<-expr_matrix[intersect(coef$gene, rownames(expr_matrix)),]
rownames(expr_matrix)[1:5]
coef1<-coef[intersect(coef$gene, rownames(expr_matrix)),]
expr_matrix<-as.data.frame(t(as.matrix(expr_matrix)))
for (ii in 1:nrow(expr_matrix)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*expr_matrix[ii,jj]
  }
  expr_matrix[ii,20]<-pcd_score
}
colnames(expr_matrix)[20]<- "PCD_scores"
seurat_obj$PCD_scores <- expr_matrix$PCD_scores
pdf(height = 8, width = 8, file = "V10T31-051-A1-PCD_scores_UMAP.pdf")
SpatialFeaturePlot(seurat_obj, features = "PCD_scores",alpha = c(0.5,1))+scale_fill_viridis(option = "A")
dev.off()




#D
p<-ggplot(data1@meta.data,aes(x=cell_type,y=PCD_scores,fill = cell_type))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+ 
  theme(legend.position = 'none')+xlab('')+ylab('PCD score')+labs(fill='')+
  scale_color_npg()
ggsave(p,file="17-celltype5-PCD_scores.pdf",width = 7,height = 4) 

p<-ggplot(data1@meta.data,aes(x=cell_type1,y=PCD_scores,fill = cell_type1))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+theme_bw()+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+ 
  theme(legend.position = 'none')+xlab('')+ylab('PCD score')+labs(fill='')+
  scale_fill_manual(values = c(pal_d3("category20")(20)[1:14],pal_locuszoom()(7)[6],pal_d3("category20")(20)[15:20]))
ggsave(p,file="17-celltype20-PCD_scores.pdf",width = 7,height = 4) 


p<-ggplot(data1@meta.data,aes(x=cell_type,y=PCD_scores,fill=status))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF6600"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", " ")),label = "p.signif")+
  theme(legend.position = 'top')+xlab('')+ylab('PCD score')+labs(fill='')
ggsave(p,file="18-celltype5-PCD_scores.pdf",width = 7,height = 4) 

p<-ggplot(data1@meta.data,aes(x=cell_type1,y=PCD_scores,fill=status))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF6600"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(symnum.args = list(cutpoints = c(0,0.001, 0.01, 0.05, 1), 
                                        symbols = c("***", "**", "*", " ")),label = "p.signif")+
  theme(legend.position = 'top')+xlab('')+ylab('PCD score')+labs(fill='')
ggsave(p,file="18-celltype20-PCD_scores.pdf",width = 7,height = 4) 




#E
data <- data1@assays$RNA@data
cluster_label <- data1@meta.data
Idents(data1)<-data1$cell_type1
markers <- FindAllMarkers(data1,pos=T)
markers_sig <- markers %>% subset(avg_log2FC > 1 & p_val_adj < 0.01)
data_1 <-data[which(rownames(data)%in% unique(markers_sig$gene)),rownames(cluster_label)]
data_1<-as.data.frame(as.matrix(data_1))
back <- data_1 %>%
  apply(1, function(x) { tapply(x, cluster_label$cell_type1 %>% as.factor(), mean) }) %>%
  t() %>% as.data.frame()
addcol <- colorRampPalette(c("#9966CC","white","#FFCC33"))(100)
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])

#GSE175457
exp<-read.table("GSE175457_expr.txt")
pdata<-read.table("GSE175457_class.txt")
data1<-exp[which(rownames(exp) %in%  rownames(back)),]
result1 <- CIBERSORT(as.matrix(back), as.matrix(data1), perm = 1000, QN = T) 
result_cell <- as.data.frame(result1[,1:21])
coef1<-coef[intersect(rownames(exp),coef$gene),]
phon<-as.data.frame(t(exp[coef1$gene,pdata$sample]))
phon$class<-pdata$outcome
for (ii in 1:nrow(phon)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*phon[ii,jj]
  }
  phon[ii,21]<-pcd_score
}
colnames(phon)[21]<-"pcd_score"
library(psych)
pp <- corr.test(phon[,21], result_cell, method = "pearson", adjust = "fdr")
cor <- pp$r 
pvalue <- pp$p 
df <- melt(cor) %>% mutate(pvalue = melt(pvalue)[, 3],
                           p_signif = symnum(pvalue, corr = FALSE, na = FALSE,  
                                             cutpoints = c(0, 0.01, 1), symbols = c("*", " "))) %>% 
  set_colnames(c("genes", "cells", "r", "p", "p_signif"))
df <-df [order(df $r),]
df$cells<-factor(df$cells,levels = df$cells)
z<-signif(df$p,2)
df$pd<-"p<0.01"
df$pd[which(df$p>0.01)]<-"p>0.01"
p<-ggplot(data = df, mapping = aes(x =cells,y=r,fill = pd)) +
  geom_bar(stat = 'identity',colour = 'gray')+ scale_fill_manual(values=c("#b8f1cc","darkgray"))+
  theme(axis.text.x = element_text(angle = 60,hjust = 1))+
  theme(panel.background=element_rect(fill="white",colour="white"), 
        axis.line=element_line(colour="black"))+
  xlab(" ")+ylab('Correlation Level')+labs(fill='')+geom_text(mapping = aes(label = z))
ggsave(p,file="GSE175457_pcd_score_cells_barplot.pdf",width = 12,height = 5)
pp <- corr.test(phon[,1:19], result_cell, method = "pearson", adjust = "fdr")
cor <- pp$r 
pvalue <- pp$p 
df1 <- melt(cor) %>% mutate(pvalue = melt(pvalue)[, 3],
                            p_signif = symnum(pvalue, corr = FALSE, na = FALSE,  
                                              cutpoints = c(0, 0.01, 1), symbols = c("*", " "))) %>% 
  set_colnames(c("genes", "cells", "r", "p", "p_signif"))
rvalue <- df1 %>% select(1, 2, 3) %>% 
  pivot_wider(names_from = "cells", values_from = r) %>% column_to_rownames(var = "genes")
pvalue <- df1 %>% select(1, 2, 5) %>% pivot_wider(names_from = "cells", values_from = p_signif) %>% 
  column_to_rownames(var = "genes")
pdf("GSE175457_pcd_gene_cells_heatmap.pdf",width = 12,height = 10)
pheatmap(rvalue[,as.character(df$cells)], scale = "none", cluster_row = T, cluster_col = F, border = NA,
         display_numbers = pvalue[,as.character(df$cells)], fontsize_number = 12, number_color = "black",
         main = " ",cellwidth = 21, cellheight = 20, color = addcol)
dev.off()





#F
Idents(data1)<-data1$cell_type
average_expression <- AverageExpression(data1, features = coef1$gene)
gene_expression_mean <- as.data.frame(average_expression$RNA[coef1$gene, ])
pdf(height = 8, width = 12, file = "19-PCD-markers-celltype-heatmap.pdf")
pheatmap(scale(as.matrix(t(gene_expression_mean))),cluster_cols = F, cluster_rows = F, 
         color =colorRampPalette(colors = c("white","#FFCC33"))(100))
dev.off()
Idents(data1)<-data1$cell_type1
average_expression <- AverageExpression(data1, features = coef1$gene)
gene_expression_mean <- as.data.frame(average_expression$RNA[coef1$gene, ])
gene_expression_mean <-gene_expression_mean [,order(colnames(gene_expression_mean))]
pdf(height = 8, width = 12, file = "19-PCD-markers-celltype1-heatmap.pdf")
pheatmap(scale(as.matrix(t(gene_expression_mean))),cluster_cols = F, cluster_rows = T, 
         color =colorRampPalette(colors = c("white","#FFCC33"))(100))
dev.off()

Idents(data1)<-data1$status
average_expression <- AverageExpression(data1, features = coef1$gene)
gene_expression_mean <- as.data.frame(average_expression$RNA[coef1$gene, ])
pdf(height = 8, width = 12, file = "19-PCD-markers-status-heatmap.pdf")
pheatmap(scale(as.matrix(t(gene_expression_mean))),cluster_cols = F, cluster_rows = F, 
         color =colorRampPalette(colors = c("white","#FFCC33"))(100))
dev.off()
Idents(data1)<-data1$status_2
average_expression <- AverageExpression(data1, features = coef1$gene)
gene_expression_mean <- as.data.frame(average_expression$RNA[coef1$gene, ])
pdf(height = 8, width = 12, file = "19-PCD-markers-status_2-heatmap.pdf")
pheatmap(scale(as.matrix(t(gene_expression_mean))),cluster_cols = F, cluster_rows = F, 
         color =colorRampPalette(colors = c("white","#FFCC33"))(100))
dev.off()
