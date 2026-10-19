library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(RColorBrewer)
library(ggplot2)
library(ggthemes)
library(ggpubr)
gene<-as.data.frame(fread("PCD-genes.txt"))
colnames(gene)
bb<-as.data.frame(table(c(gene[,1],gene[,2],gene[,3],gene[,4],gene[,5],gene[,6],gene[,7],gene[,8],gene[,9],
                          gene[,10],gene[,11],gene[,12],gene[,13],gene[,14],gene[,15],gene[,16],gene[,17],gene[,18])))
bb<-bb[-1,]
#GSE175457
gset=getGEO('GSE175457',AnnotGPL=F)
exprset <-as.data.frame(t(read.table("GSE175457_geo_sysbio_rnaseq_data.txt",header =T)))
eg = mapIds(x = org.Hs.eg.db,keys =substr(rownames(exprset),1,15),keytype = "ENSEMBL",column = "SYMBOL")
exp4<-aggregate(x=exprset,by = eg %>% list(),FUN = mean)#
rownames(exp4)<-exp4$Group.1
exp4<-exp4[,-1]
pdata4<-gset[["GSE175457_series_matrix.txt.gz"]]@phenoData@data
gene<-strsplit(pdata4$title," ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])}
pdata4$title<-gene_symbol
exp4<-exp4[,pdata4$title]
colnames(exp4)<-pdata4$geo_accession
pdata4$status<-pdata4$characteristics_ch1.2
pdata4$status[which(pdata4$status %in% c("dx: control","dx: normal","dx: other"))]<-"control"
pdata4$status[which(pdata4$status %in% c("dx: ipf"))]<-"IPF"

library(limma)
group_list <-  factor(pdata4$status)
design <- model.matrix(~0+group_list)
rownames(design) = colnames(exp4)
colnames(design) <- levels(group_list)
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
fit <- lmFit(exp4,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  
tempOutput = topTable(fit2, coef=1, n=Inf)
tT3 = na.omit(tempOutput) 
tT3$id<-rownames(tT3)
tT3$logFC<--tT3$logFC
DEGs_tt1<-tT3
DEGs_tt1<-na.omit(DEGs_tt1)
DEGs_tt1$logP<--log10(DEGs_tt1$adj.P.Val )
DEGs_tt1$group<-"not-DEGs"
DEGs_tt1$group[which((DEGs_tt1$adj.P.Val<0.01)&(DEGs_tt1$logFC>1))]="up-DEGs"
DEGs_tt1$group[which((DEGs_tt1$adj.P.Val<0.01)&(DEGs_tt1$logFC<(-1)))]="down-DEGs"

#GSEA
geneList1<-tT3$logFC
names(geneList1)<-tT3$id
geneList1<-sort(geneList1,decreasing = T)
gmt<-read.gmt("c5.go.bp.v2024.1.Hs.symbols.gmt") 
gsea@result$Description<-gsub("_"," ", tolower(gsea@result$Description)) 
result2<-gsea@result[which(gsea@result$p.adjust<0.01),]
result2$Description<-substr(result2$Description,6,100)
write.csv(result2,"GSE175457_go.bp_gsea.csv",row.names=F)
colnames(result2)
bp_frb<-tolower(c("Endothelial cell apoptotic process","Macroautophagy","Regulation of autophagy","Regulation of macroautophagy","Positive regulation of autophagy","Positive regulation of macroautophagy","Autophagy of mitochondrion","Regulation of autophagy of mitochondrion",
"Autophagosome organization","Cell killing","Leukocyte mediated cytotoxicity","Regulation of cell killing","Positive regulation of interleukin 1 beta production","Positive regulation of interleukin 1 production",
"Tumor necrosis factor superfamily cytokine production","Canonical NF-kappaB signal transduction","Cellular response to reactive oxygen species","Reactive oxygen species metabolic process","Superoxide metabolic process","Positive regulation of hydrolase activity"))
p<-ggplot(result2[which(result2$Description %in% bp_frb),],aes(x = NES, y = reorder(Description,NES), size = abs(NES),colour=p.adjust)) +
  geom_point(shape = 16) +labs(x = "NES", y = "Description")+  
  scale_colour_continuous(name="p.adjust",low="#663399",high="#FF9900")+scale_radius(range=c(2,4),name="Size")+
  guides( color = guide_colorbar(order = 1), size = guide_legend(order = 2))+theme_bw()                     
ggsave(p,file="GSE175457_go.bp_gsea_PCD.pdf",width = 6,height =6)
bp1<-setdiff(result2$Description[which(result2$NES<0)],bp_frb)
p<-ggplot(result2[which(result2$Description %in% bp1[1:20]),],aes(x = NES, y = reorder(Description,NES), size = abs(NES),colour=p.adjust)) +
  geom_point(shape = 16) +labs(x = "NES", y = "Description")+  
  scale_colour_continuous(name="p.adjust",low="#663399",high="#FF9900")+scale_radius(range=c(2,4),name="Size")+
  guides( color = guide_colorbar(order = 1), size = guide_legend(order = 2))+theme_bw()                     
ggsave(p,file="GSE175457_go.bp_gsea_downregulation.pdf",width = 5.5,height =5)
bp1<-setdiff(result2$Description[which(result2$NES>0)],bp_frb)
p<-ggplot(result2[which(result2$Description %in% bp1[1:20]),],aes(x = NES, y = reorder(Description,NES), size = abs(NES),colour=p.adjust)) +
  geom_point(shape = 16) +labs(x = "NES", y = "Description")+  
  scale_colour_continuous(name="p.adjust",low="#663399",high="#FF9900")+scale_radius(range=c(2,4),name="Size")+
  guides( color = guide_colorbar(order = 1), size = guide_legend(order = 2))+theme_bw()                     
ggsave(p,file="GSE175457_go.bp_gsea_upregulation.pdf",width = 5.3,height =5)


#heatmap
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
library(enrichplot)
library(ggridges)
library(grid)
library(pheatmap)
library(ggsci)
deg<-unique( DEGs_tt1[which(DEGs_tt1$adj.P.Val<0.01 & abs(DEGs_tt1$logFC)>1),7])
PCD_genes<-readRDS("PCD_genes.rds")
tT4<-tT3[which(tT3$id %in% PCD_genes$gene),] 
deg3<-unique(tT3[which(tT3$P.Value<0.01 & abs(tT3$logFC)>1),7])
PCD_DEG<-unique(tT4[which(tT4$P.Value<0.01 & abs(tT4$logFC)>1),7])
aa<-tT4[which(tT4$P.Value<0.01 & abs(tT4$logFC)>1),]
aa[,1:6]<-signif(aa[,1:6],digits = 3)

deg3<-unique(tT3[which(tT3$P.Value<0.01 & abs(tT3$logFC)>0.8),7])
PCD_Type<-PCD_genes[which(PCD_genes$gene %in% intersect(PCD_genes$gene,deg3)),]
unique(PCD_Type$PCD_type)
bb<-as.data.frame(table(PCD_Type$gene))
bb1<-as.character(bb$Var1[which(bb$Freq ==1)])
bb2<-as.character(bb$Var1[which(bb$Freq ==2)])
bb3<-as.character(bb$Var1[which(bb$Freq ==3)])
bb4<-as.character(bb$Var1[which(bb$Freq ==4)])
bb5<-as.character(bb$Var1[which(bb$Freq ==5)])
data2<-rbind(exp4[bb1,],exp4[bb2,],exp4[bb2,],exp4[bb3,],exp4[bb3,],exp4[bb3,],
             exp4[bb4,],exp4[bb4,],exp4[bb4,],exp4[bb4,],
             exp4[bb5,],exp4[bb5,],exp4[bb5,],exp4[bb5,],exp4[bb5,])
gene196<-data.frame(gene=c(bb1,bb2,bb2,bb3,bb3,bb3,bb4,bb4,bb4,bb4,bb5,bb5,bb5,bb5,bb5),
                    newname=paste0("gene",1:345))
rownames(data2)<-gene196$newname
gene196<-gene196[order(gene196$gene),]
PCD_Type<-PCD_Type[order(PCD_Type$gene),]
gene196$PCD_type<-PCD_Type$PCD_type
data2<-data2[gene196$newname,]
IPF<-pdata4$geo_accession[which(pdata4$status %in% "IPF")]
control<-pdata4$geo_accession[which(pdata4$status %in% "control")]
data2<-data2[,c(IPF,control)]
pdata4<-pdata4[c(IPF,control),]
annotation_col = data.frame(group =pdata4$status)
annotation_row = data.frame(PCD_Type =gene196$PCD_type)
rownames(annotation_col) = colnames(data2)
rownames(annotation_row) = rownames(data2)
unique(PCD_Type$PCD_type)
ann_colors = list(group  = c(IPF = "#FFC0CB",control = "#87CEFA"),
                  PCD_Type =c('Lysosome-dependentcelldeath'=pal_d3("category20")(20)[1+4],Autophagy=pal_d3("category20")(20)[2+4],
                              Ferroptosis=pal_d3("category20")(20)[3+4], Apoptosis=pal_d3("category20")(20)[4+4],
                              NETosis=pal_d3("category20")(20)[5+4],    Pyroptosis=pal_d3("category20")(20)[6+4],              
                              Anoikis=pal_d3("category20")(20)[7+4],    Paraptosis=pal_d3("category20")(20)[8+4],
                              Necroptosis=pal_d3("category20")(20)[9+4],Netoticcelldeath=pal_d3("category20")(20)[10+4],
                              'Immunogenic cell death'=pal_d3("category20")(20)[11+4],Cuproptosis=pal_d3("category20")(20)[12+4],
                              Methuosis=pal_d3("category20")(20)[13+4],  Entosis=pal_d3("category20")(20)[14+4],
                              Entoticcelldeath=pal_d3("category20")(20)[15+4],Parthanatos=pal_d3("category20")(20)[16+4]))
p1<-pheatmap(data2,cluster_col=F,cluster_rows = T,  show_rownames = F,show_colnames = F, 
             annotation_col = annotation_col,annotation_row = annotation_row,
             annotation_colors = ann_colors,
             colorRampPalette(c('#663399','#ffffff','#FF9900'))(100))
ggsave(p1,file="GSE175457_PCD_pheatmap.pdf",width = 6,height = 10)
dev.off()


#GSVA
library(ComplexHeatmap)
library(ggplot2)
library(corrgram)
library(ggthemes)
library(ggpubr)
library(GSEABase)
PCD_gene_type<-readRDS("PCD_gene_type.rds")
ssgsea_par <- ssgseaParam(as.matrix(exp4), PCD_gene_type)  
ssgsea_scores <- gsva(ssgsea_par)
cb2<-as.data.frame(t(ssgsea_scores))
immu<- rep(colnames(cb2),each=nrow(cb2)) #
immu <- factor(immu)
a<-pdata4$status
group <- rep(a,ncol(cb2)) 
group <- factor(group) 
value <- c() 
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value) 
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#f1f1b8","#d08a8a"))+
  theme(panel.grid = element_blank())+
  stat_compare_means(symnum.args = list(cutpoints = c(0,  0.01, 1), 
                                        symbols = c("***", " ")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('PCD Score')+labs(fill='')
ggsave(p,file="PCD_Score_barplot_GSE175457.pdf",width = 15,height = 8)


#hallmark
geneset <- getGmt(file.path("h.all.v2024.1.Hs.symbols.gmt"))
deg_exp<-as.matrix(exp4)
ssgsea_par <- ssgseaParam(deg_exp, geneset)  # all other values are default values
es.max1 <- gsva(ssgsea_par)
rownames(es.max1)<-gsub("_"," ",tolower(rownames(es.max1)))

group_list <-  factor(pdata4$status)
design <- model.matrix(~0+group_list)
rownames(design) = colnames(es.max1)
colnames(design) <- levels(group_list)
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
fit <- lmFit(es.max1,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##
fit2 <- eBayes(fit2)  ## default no trend !!!
tempOutput = topTable(fit2, coef=1, n=Inf)
tT3 = na.omit(tempOutput) 
tT3$id<-substr(tolower(rownames(tT3)),10,200)
tT3$logFC<--tT3$logFC
tT3<-tT3[order(tT3$logFC),]
tT3$id<-gsub("_"," ",tT3$id)
tT3$id<-factor(tT3$id,levels = tT3$id)

cb2<-as.data.frame(t(es.max1[rownames(tT3),]))
immu<- rep(colnames(cb2),each=nrow(cb2))
immu<-substr(tolower(immu),10,1000)
immu<-gsub("_"," ",tolower(immu))
immu <- factor(immu,levels = tT3$id) 
a<-pdata4$status
group <- rep(a,ncol(cb2))
group <- factor(group) 
value <- c() 
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#f1f1b8","#d08a8a"))+
  theme(panel.grid = element_blank())+
  stat_compare_means(symnum.args = list(cutpoints = c(0,  0.01, 1), 
                                        symbols = c("*", " ")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('GSVA Score')+labs(fill='')
ggsave(p,file="hallmark_barplot_GSE175457.pdf",width = 13,height = 8)




#immune cell
cell_gene<-read.table("Immune marker genes.txt",sep=",")
cell<-unique(as.character(cell_gene[,2]))
cell_gene_type=list()
for (i in 1:length(cell)){
  b=list(as.character(cell_gene[(as.character(cell_gene$V2))==cell[i],1]))
  cell_gene_type[i]=(b)
  names(cell_gene_type)[i]=cell[i]
}
ssgsea_par <- ssgseaParam(deg_exp, cell_gene_type)  # all other values are default values
gsva_es1 <- gsva(ssgsea_par)
group_list <-  factor(pdata4$status)
design <- model.matrix(~0+group_list)
rownames(design) = colnames(gsva_es1)
colnames(design) <- levels(group_list)
contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
fit <- lmFit(gsva_es1,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  
tempOutput = topTable(fit2, coef=1, n=Inf)
tT3 = na.omit(tempOutput) 
tT3$id<-rownames(tT3)
tT3$logFC<--tT3$logFC
tT3<-tT3[order(tT3$logFC),]
tT3$rank<-1:28
length(tT3$id[which(tT3$P.Value<0.01)])
tT3$id<-factor(tT3$id,levels = tT3$id)

cb2<-as.data.frame(t(gsva_es1[rownames(tT3),]))
immu<- rep(colnames(cb2),each=nrow(cb2))
immu <- factor(immu,levels = tT3$id) 
a<-pdata4$status
group <- rep(a,ncol(cb2)) 
group <- factor(group) 
value <- c() 
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value) 
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#f1f1b8","#d08a8a"))+
  theme(panel.grid = element_blank())+
  stat_compare_means(symnum.args = list(cutpoints = c(0,  0.01, 1), 
                                        symbols = c("*", " ")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('GSVA Score')+labs(fill='')
ggsave(p,file="Immune_barplot_GSE175457.pdf",width = 12,height = 8)

