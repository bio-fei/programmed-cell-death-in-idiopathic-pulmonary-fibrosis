model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
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
library(gg.gap)
set.seed(7889)
plan(multisession, workers = 4)
options(future.globals.maxSize = 450*1000 * 1024^2)



#A
data_seurat<-readRDS("F:\\1-PCD\\5-tissue-scRNA\\clustered_data_sub_celltype.rds")
data_seurat$status_2<-data_seurat$status_1
data_seurat$status_2[which(data_seurat$status_1 %in% "Control_LOW")]<-"Control"
data_seurat$status_2[which(data_seurat$status_1 %in% "Control_UP")]<-"Control"
Idents(data_seurat)<-data_seurat$cell_type
pdf("1-CD24-VlnPlot-major.pdf",width = 6,height = 3)
VlnPlot(data_seurat, features = c("CD24","FAIM"),pt.size = 0,ncol = 2,cols =pal_npg()(5))
dev.off()
Idents(data_seurat)<-data_seurat$cell_type1
pdf("1-CD24-VlnPlot-celltype.pdf",width = 12,height = 3)
VlnPlot(data_seurat, features = c("CD24","FAIM"),pt.size = 0,ncol = 2,cols =c(pal_d3("category20")(20)[1:14],pal_locuszoom()(7)[6],pal_d3("category20")(20)[15:20]))
dev.off()


#B
data_seurat@meta.data$cd24<-data_seurat@assays$RNA@data["CD24",]
Metadata_1<-data_seurat@meta.data
p<-ggplot(Metadata_1,aes(x=status,y=cd24,fill=status))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF9900"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(method = "t.test")+
  theme(legend.position = 'top')+xlab('')+ylab('CD24 gene expression')+labs(fill='')
ggsave(p,file="2-all-status-CD24.pdf",width = 5,height = 4) 
my_compare<-list(c('Control','IPF_UP'), c('IPF_UP','IPF_LOW'),c('Control','IPF_LOW'))
p<-ggplot(Metadata_1,aes(x=status_2,y=cd24,fill=status_2))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF9900","#CC3300"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(comparisons=my_compare,method = "t.test")+
  theme(legend.position = 'top')+xlab('')+ylab('CD24 gene expression')+labs(fill='')
ggsave(p,file="2-all-status2-CD24.pdf",width = 7,height = 4) 

Idents(data_seurat)<-data_seurat$cell_type1
#Epithelial
Epithelial_cell<-c("Goblet","Ciliated","AT2","Club","Basal" ,"AT1" )
selected_c <- WhichCells(data_seurat, idents=Epithelial_cell) 
Epithelial<- subset(data_seurat, cells = selected_c)
table(Epithelial$patient)
Epithelial@meta.data$cd24<-Epithelial@assays$RNA@data["CD24",]
Metadata_1<-Epithelial@meta.data
p<-ggplot(Metadata_1,aes(x=status,y=cd24,fill=status))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF9900"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(method = "t.test")+
  theme(legend.position = 'top')+xlab('')+ylab('CD24 gene expression')+labs(fill='')
ggsave(p,file="2-Epi-status-CD24.pdf",width = 5,height = 4) 
my_compare<-list(c('Control','IPF_UP'), c('IPF_UP','IPF_LOW'),c('Control','IPF_LOW'))
p<-ggplot(Metadata_1,aes(x=status_2,y=cd24,fill=status_2))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF9900","#CC3300"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(comparisons=my_compare,method = "t.test")+
  theme(legend.position = 'top')+xlab('')+ylab('CD24 gene expression')+labs(fill='')
ggsave(p,file="2-Epi-status2-CD24.pdf",width = 7,height = 4) 
Idents(data_seurat)<-data_seurat$cell_type1
#Ciliated
unique(data_seurat$cell_type1)
Ciliated_cell<-"Ciliated"
selected_c <- WhichCells(data_seurat, idents=Ciliated_cell) 
Ciliated<- subset(data_seurat, cells = selected_c)
table(Ciliated$patient)
Ciliated@meta.data$cd24<-Ciliated@assays$RNA@data["CD24",]
Metadata_1<-Ciliated@meta.data
p<-ggplot(Metadata_1,aes(x=status,y=cd24,fill=status))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF9900"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(method = "t.test")+
  theme(legend.position = 'top')+xlab('')+ylab('CD24 gene expression')+labs(fill='')
ggsave(p,file="2-Ciliated-status-CD24.pdf",width = 5,height = 4) 
my_compare<-list(c('Control','IPF_UP'), c('IPF_UP','IPF_LOW'),c('Control','IPF_LOW'))
p<-ggplot(Metadata_1,aes(x=status_2,y=cd24,fill=status_2))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#66CCFF","#FF9900","#CC3300"))+
  theme(panel.grid = element_blank())+
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+ 
  stat_compare_means(comparisons=my_compare,method = "t.test")+
  theme(legend.position = 'top')+xlab('')+ylab('CD24 gene expression')+labs(fill='')
ggsave(p,file="2-Ciliated-status2-CD24.pdf",width = 7,height = 4) 




#C
library(ggplot2)
logisticmodel<-readRDS("logisticmodel.rds")
coef<-data.frame(gene=names(logisticmodel[["glmBoost+Ridge"]][["coefficients"]]),
                 coef=logisticmodel[["glmBoost+Ridge"]][["coefficients"]])
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(ggplot2)
library(corrgram)
library(ggthemes)
library(ggpubr)
PCD_gene<-coef$gene[2:20]

#GSE175457
gset=getGEO('GSE175457', AnnotGPL=F)
exprset <-as.data.frame(t(read.table(" GSE175457_geo_sysbio_rnaseq_data.txt",header =T)))
eg = mapIds(x = org.Hs.eg.db,keys =substr(rownames(exprset),1,15),keytype = "ENSEMBL",column = "SYMBOL")
length(unique(eg))
exp4<-aggregate(x=exprset,by = eg %>% list(),FUN = mean) 
rownames(exp4)<-exp4$Group.1
exp4<-exp4[,-1]
pdata4<-gset[["GSE175457_series_matrix.txt.gz"]]@phenoData@data
gene<-strsplit(pdata4$title," ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])
}
length(unique(gene_symbol))
pdata4$title<-gene_symbol
exp4<-exp4[,pdata4$title]
colnames(exp4)<-pdata4$geo_accession
pdata4$status<-pdata4$characteristics_ch1.3
pdata4$status[which(pdata4$status %in% c("case_control_status: control"))]<-"0"
pdata4$status[which(pdata4$status %in% c("case_control_status: case"))]<-"1"
exp4<-exp4[PCD_gene,rownames(pdata4)]
cb2<-as.data.frame(t(exp4))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
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
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid = element_blank())+ggtitle("GSE175457")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE175457.pdf",width = 15,height = 8)


# GSE32537
gset=getGEO('GSE32537',destdir =" ",AnnotGPL=F)
exprset <- data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE32537_series_matrix.txt.gz"]]@featureData@data
gene<-strsplit(anno$gene_assignment," // ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][2])
}
exp=aggregate(x=exprset,by = gene_symbol %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE32537_series_matrix.txt.gz"]]@phenoData@data
pdata<-pdata[which(pdata$`final diagnosis:ch1` %in% c("control","IPF/UIP")),]
pdata$status<-pdata$`final diagnosis:ch1` 
pdata$status[which(pdata$status %in% c("IPF/UIP"))]<-"IPF"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid = element_blank())+ggtitle("GSE32537")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE32537.pdf",width = 15,height = 8)

#GSE47460
gset=getGEO('GSE47460',destdir =" ",AnnotGPL=F)
data1 <- data.frame(exprs(gset[[1]])) 
expMatrix <- gset[[1]] 
pdata_1<-pData(expMatrix) 
gpl1<-gset[["GSE47460-GPL14550_series_matrix.txt.gz"]]@featureData@data
gpl1<-gpl1[rownames(data1),]
exp_1<-aggregate(x=data1[,1:ncol(data1)],by =gpl1$GENE_SYMBOL  %>% list(),FUN = mean) 
rownames(exp_1)<-exp_1$Group.1
exp_1<-exp_1[-1,-1] 
a<-c(grep("IPF",c(pdata_1$characteristics_ch1.3)),
     grep("IPF",c(pdata_1$characteristics_ch1.4)),
     grep("IPF",c(pdata_1$characteristics_ch1.5)),
     grep("IPF",c(pdata_1$characteristics_ch1.6)),
     grep("IPF",c(pdata_1$characteristics_ch1.7)),
     grep("IPF",c(pdata_1$characteristics_ch1.8)),
     grep("IPF",c(pdata_1$characteristics_ch1.9)),
     grep("IPF",c(pdata_1$characteristics_ch1.10)),
     grep("IPF",c(pdata_1$characteristics_ch1.11)))
b<-c(grep("Control",c(pdata_1$characteristics_ch1)))
exp<-cbind(exp_1[,pdata_1$geo_accession[a]],exp_1[,pdata_1$geo_accession[b]])
pdata<-pdata_1
pdata<-pdata[colnames(exp),]
pdata$status<-c(rep("IPF",length(a)),rep("control",length(b)))
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid = element_blank())+ggtitle("GSE47460-GPL14550")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE47460-GPL14550.pdf",width = 15,height = 8)

data2 <- data.frame(exprs(gset[[2]])) 
expMatrix <- gset[[2]] 
pdata_2<-pData(expMatrix) 
gpl2<-gset[["GSE47460-GPL6480_series_matrix.txt.gz"]]@featureData@data
gpl2<-gpl2[rownames(data2),]
exp_2<-aggregate(x=data2[,1:ncol(data2)],by =gpl2$GENE_SYMBOL  %>% list(),FUN = mean) 
rownames(exp_2)<-exp_2$Group.1
exp_2<-exp_2[-1,-1]  #15180   153
c<-c(grep("IPF",c(pdata_2$characteristics_ch1.3)),
     grep("IPF",c(pdata_2$characteristics_ch1.4)),
     grep("IPF",c(pdata_2$characteristics_ch1.5)),
     grep("IPF",c(pdata_2$characteristics_ch1.6)),
     grep("IPF",c(pdata_2$characteristics_ch1.7)),
     grep("IPF",c(pdata_2$characteristics_ch1.8)),
     grep("IPF",c(pdata_2$characteristics_ch1.9)),
     grep("IPF",c(pdata_2$characteristics_ch1.10)),
     grep("IPF",c(pdata_2$characteristics_ch1.11)))
d<-c(grep("Control",c(pdata_2$characteristics_ch1)))
exp<-cbind(exp_2[,pdata_2$geo_accession[c]],exp_2[,pdata_2$geo_accession[d]])
pdata<-pdata_2
pdata<-pdata[colnames(exp),]
pdata$status<-c(rep("IPF",length(c)),rep("control",length(d)))
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid = element_blank())+ggtitle("GSE47460-GPL6480")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE47460-GPL6480.pdf",width = 15,height = 8)

 
#GSE150910
gset=getGEO('GSE150910',destdir =" ",AnnotGPL=F)
exprset <-fread(" GSE150910_gene-level_count_file.csv")
rownames(exprset)<-exprset$symbol
exp<-as.data.frame(log2(exprset[,-1]+1))
rownames(exp)<-exprset$symbol
pdata<-gset[["GSE150910_series_matrix.txt.gz"]]@phenoData@data
pdata<-pdata[which(pdata$characteristics_ch1.3 %in% c("diagnosis: control","diagnosis: ipf")),]
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status<-pdata$characteristics_ch1.3
pdata$status[which(pdata$status %in% c("diagnosis: ipf"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diagnosis: control"))]<-"control"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE150910")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE150910.pdf",width = 15,height = 8)


 
#GSE134692
gset=getGEO('GSE134692',destdir =" ",AnnotGPL=F)
exprset <-read.table(" GSE134692_tmm_normalized_and_filtered_log2CPM.txt",header = T)
eg = mapIds(x = org.Hs.eg.db,keys =exprset$ensembl,keytype = "ENSEMBL",column = "SYMBOL")
exp<-aggregate(x=exprset[,2:ncol(exprset)],by = eg %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE134692_series_matrix.txt.gz"]]@phenoData@data
table(pdata$characteristics_ch1.2)
pdata<-pdata[which(pdata$characteristics_ch1.2 %in% c("diseasestatus: IPF","diseasestatus: Normal")),]
pdata$status<-pdata$characteristics_ch1.2 
pdata$status[which(pdata$status %in% c("diseasestatus: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diseasestatus: Normal"))]<-"control"
exp<-exp[,pdata$description]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE134692")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE134692.pdf",width = 15,height = 8)

#GSE124685
gset=getGEO('GSE124685',destdir =" ",AnnotGPL=F)
exprset <-read.table(" GSE124685_FPKM_Matrix_84_samples.txt",header = T)
exp<-aggregate(x=exprset[,9:ncol(exprset)],by = exprset$gene_id %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-log2(exp[-1,-1]+1)
pdata<-gset[["GSE124685_series_matrix.txt.gz"]]@phenoData@data
table(pdata$characteristics_ch1)
pdata$status<-pdata$characteristics_ch1
pdata$status[which(pdata$status %in% c("disease state: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("disease state: healthy donor"))]<-"control"
exp<-exp[,pdata$geo_accession]
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE124685")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE124685.pdf",width = 15,height = 8)

 
#GSE110147
gset=getGEO('GSE110147',destdir =" ",AnnotGPL=F)
exprset <- data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE110147_series_matrix.txt.gz"]]@featureData@data
gene<-strsplit(anno$gene_assignment," // ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][2])
}
exp=aggregate(x=exprset,by = gene_symbol %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[,-1]
pdata<-gset[["GSE110147_series_matrix.txt.gz"]]@phenoData@data
table(pdata$characteristics_ch1)
pdata<-pdata[which(pdata$characteristics_ch1 %in% c("disease state: Idiopathic pulmonary fibrosis",
                                                    "disease state: Normal control")),]
pdata$status<-pdata$characteristics_ch1
pdata$status[which(pdata$status %in% c("disease state: Idiopathic pulmonary fibrosis"))]<-"IPF"
pdata$status[which(pdata$status %in% c("disease state: Normal control"))]<-"control"
exp<-exp[,pdata$geo_accession]
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE110147")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE110147.pdf",width = 15,height = 8)

 
#GSE53845
gset=getGEO('GSE53845',destdir =" ",AnnotGPL=F)
exprset <- data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE53845_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE53845_series_matrix.txt.gz"]]@phenoData@data
table(pdata$characteristics_ch1.2)
pdata$status<-pdata$characteristics_ch1.2
pdata$status[which(pdata$status %in% c("diagnosis: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diagnosis: control"))]<-"control"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE53845")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE53845.pdf",width = 15,height = 8)

 
#GSE92592
gset=getGEO('GSE92592',destdir =" ",AnnotGPL=F)
exprset <- read.table(" GSE92592_gene.counts.txt",header =T)
gene<-strsplit(rownames(exprset),".c")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])
}
length(unique(gene_symbol))
exp=aggregate(x=exprset,by = gene_symbol %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-log2(exp[,-1]+1)
pdata<-gset[["GSE92592_series_matrix.txt.gz"]]@phenoData@data
table(pdata$`disease state:ch1`)
pdata$status<-pdata$`disease state:ch1`
pdata$status[which(pdata$status %in% c("Idiopathic Pulmonary Fibrosis"))]<-"IPF"
pdata$status[which(pdata$status %in% c("Control"))]<-"control"
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE92592")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE92592.pdf",width = 15,height = 8)
 
#GSE231693
gset=getGEO('GSE231693',destdir =" ",AnnotGPL=F)
exprset <- read.table(" GSE231693_count.tsv",header =T)
eg = mapIds(x = org.Hs.eg.db,keys =rownames(exprset),keytype = "ENSEMBL",column = "SYMBOL")
length(unique(eg))
exp<-aggregate(x=exprset,by = eg %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-log2(exp[,-1]+1)
pdata<-gset[["GSE231693_series_matrix.txt.gz"]]@phenoData@data
table(pdata$characteristics_ch1.1)
pdata<-pdata[which(pdata$characteristics_ch1.1 %in% c("disease state: idiopathic pulmonary fibrosis","disease state: normal")),]
pdata$status<-pdata$characteristics_ch1.1
pdata$status[which(pdata$status %in% c("disease state: idiopathic pulmonary fibrosis"))]<-"IPF"
pdata$status[which(pdata$status %in% c("disease state: normal"))]<-"control"
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE231693")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE231693.pdf",width = 15,height = 8)

 
#GSE213001
gset=getGEO('GSE213001',destdir =" ",AnnotGPL=F)
exprset <- as.data.frame(fread(" GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv",header =T))
eg = mapIds(x = org.Hs.eg.db,keys =exprset$V1,keytype = "ENSEMBL",column = "SYMBOL")
length(unique(eg))
exp<-aggregate(x=exprset[,-1],by = eg %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-log2(exp[,-1]+1)
pdata<-gset[["GSE213001_series_matrix.txt.gz"]]@phenoData@data
table(pdata$characteristics_ch1.10)
pdata<-pdata[which(pdata$characteristics_ch1.10 %in% c("diagnosis: IPF","diagnosis: NDC")),]
pdata$status<-pdata$characteristics_ch1.10
pdata$status[which(pdata$status %in% c("diagnosis: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diagnosis: NDC"))]<-"control"
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE213001")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE213001.pdf",width = 15,height = 8)

 
#GSE173355
gset=getGEO('GSE173355',destdir =" ",AnnotGPL=F)
exprset <-as.data.frame(read.table(" GSE173355_PolyA_Counts.txt",header =T))
gene<-strsplit(exprset$target_id,"\\|")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][2])
}
length(unique(gene_symbol))
eg = mapIds(x = org.Hs.eg.db,keys =substr(gene_symbol,1,15),keytype = "ENSEMBL",column = "SYMBOL")
length(unique(eg))
exp<-aggregate(x=exprset[,-1],by = eg %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-log2(exp[,-1]+1)
pdata<-gset[["GSE173355_series_matrix.txt.gz"]]@phenoData@data
table(pdata$characteristics_ch1.4)
pdata$status<-pdata$characteristics_ch1.4
pdata$status[which(pdata$status %in% c("diagnosis: ipf"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diagnosis: normal"))]<-"control"
pdata$title<-substr(pdata$title,16,100)
pdata$title<-gsub("-","\\.",pdata$title)
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE173355")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE173355.pdf",width = 15,height = 8)

 
#GSE99621
gset=getGEO('GSE99621',destdir =" ",AnnotGPL=F)
exprset <-read.table(" GSE99621_datamatrix.txt",header =T)
rownames(exprset)<-exprset$Symbol
exp<-log2(exprset[,-1]+1)
pdata<-gset[["GSE99621_series_matrix.txt.gz"]]@phenoData@data
table(pdata$characteristics_ch1.1)
pdata$status<-pdata$characteristics_ch1.1
pdata$status[which(pdata$status %in% c("disease_status: idiopathic pulmonary fibrosis"))]<-"IPF"
pdata$status[which(pdata$status %in% c("disease_status: healthy control"))]<-"control"
exp<-exp[,pdata$description.1]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE99621GSE99621")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE99621.pdf",width = 15,height = 8)

 
#GSE199949
gset=getGEO('GSE199949',destdir =" ",AnnotGPL=F)
exprset <-as.data.frame(fread(" GSE199949_normalized.txt",header =T))
rownames(exprset)<-exprset$SymboP
exprset <-exprset[,10:51]
exp<-as.data.frame(log2(exprset+1))
pdata<-gset[["GSE199949_series_matrix.txt.gz"]]@phenoData@data
table(pdata$`diagnosis:ch1`)
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status<-pdata$`diagnosis:ch1`
pdata$status[which(pdata$status %in% c("Donor"))]<-"control"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE199949")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE199949.pdf",width = 15,height = 8)


#GSE70866
gset=getGEO('GSE70866',destdir =" ",AnnotGPL=F)
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE70866-GPL14550_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE70866-GPL14550_series_matrix.txt.gz"]]@phenoData@data
table(pdata$`diagnosis:ch1`)
exp<-exp[,pdata$geo_accession]
pdata$status<-pdata$`diagnosis:ch1`
pdata$status[which(pdata$status %in% c("healthy donor"))]<-"control"
#pdata$status[which(pdata$status %in% c("diagnosis: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE70866")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE70866.pdf",width = 15,height = 8)


#GSE33566
gset=getGEO('GSE33566',destdir =" ",AnnotGPL=F) 
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE33566_series_matrix.txt.gz"]]@featureData@data
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE33566_series_matrix.txt.gz"]]@phenoData@data
table(pdata$source_name_ch1)
exp<-exp[,pdata$geo_accession]
pdata$status<-pdata$source_name_ch1
pdata$status[which(pdata$status %in% c("healthy control (N)_PBB"))]<-"control"
pdata$status[which(pdata$status %in% c("idiopathic pulmonary fibrosis (IPF)_PBB"))]<-"IPF"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene,rownames(pdata)]
cb2<-as.data.frame(t(exp))
immu<- rep(colnames(cb2),each=nrow(cb2)) 
immu <- factor(immu)  
a<-pdata$status
group <- rep(a,ncol(cb2))  
group <- factor(group)  
value <- c()  
for (j in 1:ncol(cb2)) { value<-c(value,cb2[,j])}
value<-as.numeric(value)
Data <- data.frame(immu_cell=immu,group=group,value=value)  
p<-ggplot(Data,aes(x=immu_cell,y=value,fill=group))+
  geom_boxplot(width=0.7,size=0.3,outlier.color = NA)+
  theme_bw()+scale_fill_manual(values = c("#FFCC33","#9966CC"))+
  theme(panel.grid =  element_blank())+ggtitle("GSE33566")+
  stat_compare_means(symnum.args = list(cutpoints = c(0, 0.001, 0.01, 1), 
                                        symbols = c("***", "*", "ns")),label = "p.signif")+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  theme(legend.position = 'top')+xlab('')+ylab('gene expression')+labs(fill='')
ggsave(p,file="PCD_feature_exp_barplot_GSE33566.pdf",width = 15,height = 8)



#D-E
library(Seurat)
library(clustree)
library(ggsci)
file_path<-".\\E-MTAB-14121\\"
files<-list.files(file_path)
for (i in files[c(1:28)]) {
  print(i)
  data.dir<-paste0(file_path,i)
  seurat_obj <- Load10X_Spatial(data.dir =data.dir)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, "^MT-")
  seurat_obj[["percent.ribo"]] <-PercentageFeatureSet(seurat_obj, "^RP[SL]")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, "^HB[AB]")
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj,vars.to.regress = "percent.mt")   
  seurat_obj <- RunPCA(seurat_obj)
  pct <-seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)  #6
  pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:pc.use)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pc.use)
  seurat_obj <- FindClusters(object = seurat_obj,resolution = 0.5)
  pdf(height = 8, width = 18, file = paste0(i,"-EPCAM-FOXJ1-CD24_UMAP.pdf"))
  p<-SpatialFeaturePlot(seurat_obj, features = c("EPCAM","FOXJ1","CD24"),alpha = c(0.5,1))
  print(p)
  dev.off()
  mycolor <- c('white', '#ac5e74')
  pdf(height =5, width = 8, file = paste0(i,"_markers_dot.pdf"))
  p<-DotPlot(object = seurat_obj,features = c("EPCAM","FOXJ1","CD24"),cols= mycolor) + RotatedAxis()+coord_flip()
  print(p)
  dev.off()
  
}




#F-G
library(Seurat)
library(clustree)
library(ggsci)
file_path<-".\\GSE267904_RAW\\"
files<-list.files(file_path)
for (i in files[1:24]) {
  print(i)
  data.dir<-paste0(file_path,i)
  seurat_obj <- Read10X(data.dir =data.dir)
  image2 <- Read10X_Image(image.dir = file.path(data.dir, "spatial"), filter.matrix = TRUE)
  seurat_obj <- CreateSeuratObject(counts = seurat_obj, assay = "Spatial")
  image2 <- image2[Cells(x = seurat_obj)]
  DefaultAssay(seurat_obj = image2) <- "Spatial"
  seurat_obj[["slice1"]] <- image2
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  seurat_obj[["percent.ribo"]] <-PercentageFeatureSet(seurat_obj, "^Rp[sl]")
  seurat_obj[["percent.hb"]] <- PercentageFeatureSet(seurat_obj, "^Hb[ab]")
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj,vars.to.regress = "percent.mt")    
  seurat_obj <- RunPCA(seurat_obj)
  pct <-seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  cumu <- cumsum(pct)  #6
  pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:pc.use)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:pc.use)
  seurat_obj <- FindClusters(object = seurat_obj,resolution = 0.5)
  pdf(height = 12, width = 12, file = paste0(i,"-EPCAM-FOXJ1-CD24_UMAP.pdf"))
  p<-SpatialFeaturePlot(seurat_obj, features = c("Krt8","Krt18","Cd24a","Lrrc23","Foxj1"),alpha = c(0.5,1))
  print(p)
  dev.off()
  mycolor <- c('white', '#ac5e74')
  pdf(height =5, width = 8, file = paste0(i,"_markers_dot.pdf"))
  p<-DotPlot(object = seurat_obj,features = c("Krt8","Krt18","Cd24a","Lrrc23","Foxj1"),cols= mycolor) + RotatedAxis()+coord_flip()
  print(p)
  dev.off()
  
}



#H
#GSE175457
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
data<-read.table("Training_expr.txt")
data_calss<-read.table("Training_class.txt",header = T)
phon<-as.data.frame(t(data[coef$gene,order(colnames(data))]))
data_calss<-data_calss[order(data_calss$sample),]
phon$class<-data_calss$outcome
coef1<-coef[intersect(rownames(data),coef$gene),]
for (ii in 1:nrow(phon)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*phon[ii,jj]
  }
  phon[ii,21]<-pcd_score
}
colnames(phon)[21]<-"pcd_score"
library(ggplot2)
library(ggthemes)
library(ggpubr)
p <- ggplot(phon, aes(CD24,pcd_score)) +
  geom_point(aes(col=as.character(class))) +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#f1f1b8","#d08a8a"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=phon, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="PCD_score-CD24.pdf",width =5,height =3)



#I
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(survival)
library(survminer)
library(forestplot)
#GSE93606
gset=getGEO('GSE93606', AnnotGPL=F) 
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE93606_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
gene<-strsplit(anno$gene_assignment," // ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][2])
}
exp=aggregate(x=exprset,by = gene_symbol %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[,-1]
exp<-exp[intersect(rownames(exp),coef$gene),]
pdata<-gset[["GSE93606_series_matrix.txt.gz"]]@phenoData@data
pdata<-pdata[which(pdata$`disease state:ch1` %in% "Idiopathic Pulmonary Fibrosis"),]
surv<-as.data.frame(t(exp[,rownames(pdata)]))
surv$status<-as.numeric(pdata$`composite_end_point:ch1`)
surv$time<-as.numeric(pdata$`survival (months):ch1`)
surv<-surv[complete.cases(surv),]
res.cut <- surv_cutpoint(surv, time = "time", event = "status", variables ="CD24")
res.cat <- surv_categorize(res.cut)
surv$group<-res.cat[,3]
fit <- survfit(Surv(time,status) ~ group, data =surv)
p<-ggsurvplot(fit, data =surv, pval = T,  conf.int = F,risk.table = F, 
              legend.title="CD24", title="GSE93606 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE93606_CD24_os_surv.pdf",width = 9,height =6.6)
print(p,newpage = F)
dev.off()

#GSE28042
gset=getGEO('GSE28042', AnnotGPL=F) 
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE28042_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
dim(exp)    #18614    94  75/19
exp<-exp[intersect(rownames(exp),coef$gene),]
pdata<-gset[["GSE28042_series_matrix.txt.gz"]]@phenoData@data
pdata<-pdata[which(pdata$`disease status` %in% "idiopathic pulmonary fibrosis (IPF)"),]
gene<-strsplit(pdata$`time to outcome:ch1`," ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])
}
pdata$time<-gene_symbol
surv<-as.data.frame(t(exp[,rownames(pdata)]))
surv$status<-as.numeric(pdata$`outcome:ch1`)
surv$time<-as.numeric(pdata$time)
surv<-surv[complete.cases(surv),]
res.cut <- surv_cutpoint(surv, time = "time", event = "status", variables ="CD24")
res.cat <- surv_categorize(res.cut)
surv$group<-res.cat[,3]
surv$group<-"low"
surv$group[which(surv$CD24 > mean(surv$CD24) )]<-"high"
fit <- survfit(Surv(time,status) ~ group, data =surv)
p<-ggsurvplot(fit, data =surv, pval = T,  conf.int = F,risk.table = F, 
              legend.title="CD24", title="GSE28042 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE28042_CD24_os_surv.pdf",width = 9,height =6.6)
print(p,newpage = F)
dev.off()

#GSE27957
gset=getGEO('GSE27957', AnnotGPL=F) 
library(oligo)
celfiles <- list.files("GSE27957_RAW")
data.raw <- read.celfiles(filenames = file.path(celfiles))
data.raw <- rma(data.raw)　 
exprset <- as.data.frame(exprs(data.raw))
anno <-as.data.frame(fread("GSE27957_GPL5175.txt"))
anno <-anno[which(anno$ID %in% rownames(exprset)),]
gene=anno$gene_assignment %>% strsplit(" // ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][2])
}
anno$`Gene Symbol`<-gene_symbol
anno<-anno[order(anno$ID),]
exprset<-exprset[order(rownames(exprset)),]
colnames(exprset)<-substr(colnames(exprset),1,9)
exp<-aggregate(x=exprset,by = anno$`Gene Symbol` %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[,-1]
exp<-exp[intersect(rownames(exp),coef$gene),]
pdata<-fread("F:\\1-PCD\\data\\GSE27957_clinical.txt")  
surv<-as.data.frame(t(exp[,pdata$accession]))
surv$status<-as.numeric(pdata$outcome)
surv$time<-as.numeric(substr(pdata$`time to outcome`,1,5))*365
surv<-surv[complete.cases(surv),]
res.cut <- surv_cutpoint(surv, time = "time", event = "status", variables ="CD24")
res.cat <- surv_categorize(res.cut)
surv$group<-res.cat[,3]
fit <- survfit(Surv(time,status) ~ group, data =surv)
p<-ggsurvplot(fit, data =surv, pval = T,  conf.int = F,risk.table = F, 
              legend.title="CD24", title="GSE27957 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE27957_CD24_os_surv.pdf",width = 9,height =6.6)
print(p,newpage = F)
dev.off()


#GSE70866
gset=getGEO('GSE70866', AnnotGPL=F)
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE70866-GPL14550_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE70866-GPL14550_series_matrix.txt.gz"]]@phenoData@data
table(pdata$`diagnosis:ch1`)
pdata<-pdata[which(pdata$`diagnosis:ch1` %in% "IPF"),]
exp<-exp[intersect(rownames(exp),coef$gene),]
surv<-as.data.frame(t(exp[,pdata$geo_accession]))
surv$status<-as.numeric(pdata$`survival status, 0 = censored, 1 = death:ch1`)
surv$time<-as.numeric(pdata$`time to death (days):ch1`)
surv<-surv[complete.cases(surv),]
res.cut <- surv_cutpoint(surv, time = "time", event = "status", variables ="CD24")
res.cat <- surv_categorize(res.cut)
surv$group<-res.cat[,3]
fit <- survfit(Surv(time,status) ~ group, data =surv)
p<-ggsurvplot(fit, data =surv, pval = T,  conf.int = F,risk.table = F, 
              legend.title="CD24", title="GSE70866-GPL14550 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE70866-GPL14550_CD24_os_surv.pdf",width = 9,height =6.6)
print(p,newpage = F)
dev.off()


exprset <-data.frame(exprs(gset[[2]])) 
anno<-gset[["GSE70866-GPL17077_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE70866-GPL17077_series_matrix.txt.gz"]]@phenoData@data
table(pdata$`diagnosis:ch1`)
exp<-exp[intersect(rownames(exp),coef$gene),]
surv<-as.data.frame(t(exp[,pdata$geo_accession]))
surv$status<-as.numeric(pdata$`survival status, 0 = censored, 1 = death:ch1`)
surv$time<-as.numeric(pdata$`time to death (days):ch1`)
surv<-surv[complete.cases(surv),]
res.cut <- surv_cutpoint(surv, time = "time", event = "status", variables ="CD24")
res.cat <- surv_categorize(res.cut)
surv$group<-res.cat[,3]
fit <- survfit(Surv(time,status) ~ group, data =surv)
p<-ggsurvplot(fit, data =surv, pval = T,  conf.int = F,risk.table = F, 
              legend.title="CD24", title="GSE70866-GPL17077 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE70866-GPL17077_CD24_os_surv.pdf",width = 9,height =6.6)
print(p,newpage = F)
dev.off()



#J-K
ppi<-read.table("string_interactions_short.tsv",header =F,sep = "\t")
ppi<-ppi[,1:2]
colnames(ppi)<-c("gene1","gene2")
gene<-unique(c(ppi$gene1,ppi$gene2))

library(clusterProfiler)
library(org.Hs.eg.db)
library(xlsx)
library(ggplot2)
library(GOplot)
library(stringr)
library(enrichplot)
eg = mapIds(x = org.Hs.eg.db,keys = gene,keytype = "SYMBOL",column = "ENTREZID")
library(ggsci)
library(viridis)
library(RColorBrewer)
R.utils::setOption( "clusterProfiler.download.method",'auto' )
kegg <- enrichKEGG(gene    = eg, organism     = 'hsa', pvalueCutoff = 0.05)
kegg <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
kegg_result<-kegg@result
kegg_result$Description[1:20]
pdf(height =6, width = 8, file = "CD24_gene_kegg.pdf")
treemap(kegg_result[1:20,],index = "Description", 
        vSize = "Count", type = "manual",
        vColor = "p.adjust" ,palette = viridis(20))
dev.off()

