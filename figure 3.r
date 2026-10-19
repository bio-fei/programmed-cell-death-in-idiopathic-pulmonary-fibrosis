#A
logisticmodel<-readRDS("logisticmodel.rds")
coef<-data.frame(gene=names(logisticmodel[["glmBoost+Ridge"]][["coefficients"]]),
                 coef=logisticmodel[["glmBoost+Ridge"]][["coefficients"]])
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(survival)
library(survminer)


#GSE27957
library(oligo)
celfiles <- list.files("GSE27957_RAW")
setwd("GSE27957_RAW")
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
pdata<-fread("GSE27957_clinical.txt")  
pdata$status<-as.numeric(pdata$outcome)
pdata$time<-as.numeric(substr(pdata$`time to outcome`,1,5))*365
pdata<-pdata[,c(1,13,14)]
colnames(pdata)<-c("sample","OS","OS.time")
pdata<-pdata[complete.cases(pdata),]
surv<-cbind(as.data.frame(t(exp[,pdata$sample])),pdata)
library(survival)
library(survminer)
library(glmnet) 
set.seed(222)

surv<-surv[which(surv$OS.time>0),]
surv1<-survival::Surv(surv$OS.time,surv$OS)
fit<-glmnet(as.matrix(surv[,1:19]),surv1,family='cox',alpha = 1,nfolds =1000)
cvfit<-cv.glmnet(as.matrix(surv[,1:19]),surv1,family='cox',alpha = 1,nfolds =1000)
coef_lasso = coef(fit, s = cvfit$lambda.min) 
index = which(coef_lasso != 0) 
geneCoef =as.data.frame (cbind(Gene=row.names(coef_lasso)[index] ,Coef=coef_lasso[index]) )
lasso_exp<-surv[,c(which(colnames(surv)%in% geneCoef[,1]),21:22)]
for (i in 1:nrow(lasso_exp)) {
  a<-0
  for (j in 1:4) {
    a<-a+as.numeric(geneCoef[j,2])*lasso_exp[i,j]
  }
  lasso_exp[i,7]<-a
}
rownames(lasso_exp)<-rownames(surv)
res.cut <- surv_cutpoint(lasso_exp, time = "OS.time", event = "OS", variables ="V7")
res.cat <- surv_categorize(res.cut)
lasso_exp$group<-res.cat[,3]
fit <- survfit(Surv(OS.time,OS) ~ group, data =lasso_exp)
p<-ggsurvplot(fit, data =lasso_exp, pval = T,  conf.int = F,risk.table = F, 
              legend.title="pcd_score", title="GSE27957 Overall survival",
              ylab="percentage",xlab = " Time (Months)",
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE27957_os_surv.pdf",width = 9,height =6.6)
print(p,newpage = F)
dev.off()



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
pdata$status<-as.numeric(pdata$`composite_end_point:ch1`)
pdata$time<-as.numeric(pdata$`survival (months):ch1`)*30
pdata<-pdata[complete.cases(pdata),]
pdata<-pdata[,c(2,3,44)]
colnames(pdata)<-c("sample","OS","OS.time")
exp<-exp[,rownames(pdata)]
surv<-cbind(as.data.frame(t(exp[,rownames(pdata)])),pdata)
lasso_exp<-surv[,c(which(colnames(surv)%in% geneCoef[,1]),21:22)]
for (i in 1:nrow(lasso_exp)) {
  a<-0
  for (j in 1:4) {
    a<-a+as.numeric(geneCoef[j,2])*lasso_exp[i,j]
  }
  lasso_exp[i,7]<-a
}
rownames(lasso_exp)<-rownames(surv)
res.cut <- surv_cutpoint(lasso_exp, time = "OS.time", event = "OS", variables ="V7")
res.cat <- surv_categorize(res.cut)
lasso_exp$group<-res.cat[,3]
fit <- survfit(Surv(OS.time,OS) ~ group, data =lasso_exp)
p<-ggsurvplot(fit, data =lasso_exp, pval = T,  conf.int = F,risk.table = F, 
              legend.title="pcd_score", title="GSE93606 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE93606_os_surv.pdf",width = 9,height =6.6)
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
exp<-exp[intersect(rownames(exp),coef$gene),]
pdata<-gset[["GSE28042_series_matrix.txt.gz"]]@phenoData@data
pdata<-pdata[which(pdata$`disease status` %in% "idiopathic pulmonary fibrosis (IPF)"),]
gene<-strsplit(pdata$`time to outcome:ch1`," ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])
}
pdata$time<-gene_symbol
pdata$status<-as.numeric(pdata$`outcome:ch1`)
pdata$time<-as.numeric(pdata$time)*365
pdata<-pdata[,c(2,3,49)]
pdata<-pdata[complete.cases(pdata),]
exp<-exp[,pdata$geo_accession]
colnames(pdata)<-c("sample","OS","OS.time")
surv<-cbind(as.data.frame(t(exp[,rownames(pdata)])),pdata)
lasso_exp<-surv[,c(which(colnames(surv)%in% geneCoef[,1]),21:22)]
for (i in 1:nrow(lasso_exp)) {
  a<-0
  for (j in 1:4) {
    a<-a+as.numeric(geneCoef[j,2])*lasso_exp[i,j]
  }
  lasso_exp[i,7]<--a
}
rownames(lasso_exp)<-rownames(surv)
res.cut <- surv_cutpoint(lasso_exp, time = "OS.time", event = "OS", variables ="V7")
res.cat <- surv_categorize(res.cut)
lasso_exp$group<-res.cat[,3]
fit <- survfit(Surv(OS.time,OS) ~ group, data =lasso_exp)
p<-ggsurvplot(fit, data =lasso_exp, pval = T,  conf.int = F,risk.table = F, 
              legend.title="pcd_score", title="GSE28042 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE28042_os_surv.pdf",width = 9,height =6.6)
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
pdata$status<-as.numeric(pdata$`survival status, 0 = censored, 1 = death:ch1`)
pdata$time<-as.numeric(pdata$`time to death (days):ch1`)
pdata<-pdata[,c(2,3,47)]
pdata<-pdata[complete.cases(pdata),]
exp<-exp[,pdata$geo_accession]
colnames(pdata)<-c("sample","OS","OS.time")
surv<-cbind(as.data.frame(t(exp[,rownames(pdata)])),pdata)
lasso_exp<-surv[,c(which(colnames(surv)%in% geneCoef[,1]),21:22)]
for (i in 1:nrow(lasso_exp)) {
  a<-0
  for (j in 1:4) {
    a<-a+as.numeric(geneCoef[j,2])*lasso_exp[i,j]
  }
  lasso_exp[i,7]<--a
}
rownames(lasso_exp)<-rownames(surv)
res.cut <- surv_cutpoint(lasso_exp, time = "OS.time", event = "OS", variables ="V7")
res.cat <- surv_categorize(res.cut)
lasso_exp$group<-res.cat[,3]
fit <- survfit(Surv(OS.time,OS) ~ group, data =lasso_exp)
p<-ggsurvplot(fit, data =lasso_exp, pval = T,  conf.int = F,risk.table = F, 
              legend.title="pcd_score", title="GSE70866-GPL14550 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE70866-GPL14550_os_surv.pdf",width = 9,height =6.6)
print(p,newpage = F)
dev.off()


#GSE70866-GPL17077
exprset <-data.frame(exprs(gset[[2]])) 
anno<-gset[["GSE70866-GPL17077_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE70866-GPL17077_series_matrix.txt.gz"]]@phenoData@data
exp<-exp[intersect(rownames(exp),coef$gene),]
pdata$status<-as.numeric(pdata$`survival status, 0 = censored, 1 = death:ch1`)
pdata$time<-as.numeric(pdata$`time to death (days):ch1`)
pdata<-pdata[,c(2,3,47)]
pdata<-pdata[complete.cases(pdata),]
exp<-exp[,pdata$geo_accession]
colnames(pdata)<-c("sample","OS","OS.time")
surv<-cbind(as.data.frame(t(exp[,rownames(pdata)])),pdata)
lasso_exp<-surv[,c(which(colnames(surv)%in% geneCoef[,1]),21:22)]
for (i in 1:nrow(lasso_exp)) {
  a<-0
  for (j in 1:4) {
    a<-a+as.numeric(geneCoef[j,2])*lasso_exp[i,j]
  }
  lasso_exp[i,7]<-a
}
rownames(lasso_exp)<-rownames(surv)
res.cut <- surv_cutpoint(lasso_exp, time = "OS.time", event = "OS", variables ="V7")
res.cat <- surv_categorize(res.cut)
lasso_exp$group<-res.cat[,3]
fit <- survfit(Surv(OS.time,OS) ~ group, data =lasso_exp)
p<-ggsurvplot(fit, data =lasso_exp, pval = T,  conf.int = F,risk.table = F, 
              legend.title="pcd_score", title="GSE70866-GPL17077 Overall survival",
              ylab="percentage",xlab = " Time (Months)", 
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "GSE70866-GPL17077_os_surv.pdf",width = 9,height =6.6)
print(p,newpage = F)
dev.off()


#B
rm(list = ls())
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(survival)
library(survminer)

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
exp1<-exp[,-1]
pdata<-gset[["GSE93606_series_matrix.txt.gz"]]@phenoData@data
pdata1<-pdata[which(pdata$`disease state:ch1` %in% "Idiopathic Pulmonary Fibrosis"),]


#GSE28042
gset=getGEO('GSE28042', AnnotGPL=F) 
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE28042_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp2<-exp[-1,-1]
pdata<-gset[["GSE28042_series_matrix.txt.gz"]]@phenoData@data
pdata2<-pdata[which(pdata$`disease status` %in% "idiopathic pulmonary fibrosis (IPF)"),]
gene<-strsplit(pdata2$`time to outcome:ch1`," ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])
}
pdata2$time<-gene_symbol


#GSE27957
gset=getGEO('GSE27957', AnnotGPL=F) 
library(oligo)
celfiles <- list.files("GSE27957_RAW")
setwd("GSE27957_RAW")
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
exp3<-exp[,-1]
pdata3<-fread("GSE27957_clinical.txt")  


#GSE70866
gset=getGEO('GSE70866', AnnotGPL=F)
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE70866-GPL14550_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp4<-exp[-1,-1]
pdata<-gset[["GSE70866-GPL14550_series_matrix.txt.gz"]]@phenoData@data
table(pdata$`diagnosis:ch1`)
pdata4<-pdata[which(pdata$`diagnosis:ch1` %in% "IPF"),]

#####
exprset <-data.frame(exprs(gset[[2]])) 
anno<-gset[["GSE70866-GPL17077_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp5<-exp[-1,-1]
pdata5<-gset[["GSE70866-GPL17077_series_matrix.txt.gz"]]@phenoData@data
table(pdata5$`diagnosis:ch1`)

library(sva)
gene<-Reduce(intersect,list(rownames(exp1),rownames(exp2),rownames(exp3),rownames(exp4),rownames(exp5)))
exp1<-exp1[gene,pdata1$geo_accession]
exp2<-exp2[gene,pdata2$geo_accession]
exp3<-exp3[gene,pdata3$accession]
exp4<-exp4[gene,pdata4$geo_accession]
exp5<-exp5[gene,pdata5$geo_accession]
count_matrix<-(cbind(exp1,exp2,exp3,exp4,exp5))

library(FactoMineR)
library(factoextra)
count_matrix<-as.matrix(t(count_matrix))
cla <-c(rep("GSE93606",154),rep("GSE28042",75),rep("GSE27957",45),rep("GSE70866_GPL14550",112),rep("GSE70866_GPL17077",64))
iris.pca <- PCA(count_matrix,graph = F)
pdf("pca1.pdf",width = 6,height = 4.5)
ind.p<-fviz_pca_ind(iris.pca,geom.ind = "point", col.ind = cla, 
                    palette = c( "#FF9999","#FFCC66","#99CCFF","#d9b8f1","#ddff95"),
                    addEllipses = TRUE, 
                    legend.title = "" )
ggpubr::ggpar(ind.p,title = "")
dev.off()
count_matrix1<-as.matrix(t(adjusted_counts))
iris.pca <- PCA(count_matrix1,graph =F)
pdf("pca2.pdf",width = 6,height = 4.5)
ind.p<-fviz_pca_ind(iris.pca,geom.ind = "point", col.ind = cla, 
                    palette = c( "#FF9999","#FFCC66","#99CCFF","#d9b8f1","#ddff95"),
                    addEllipses = TRUE, 
                    legend.title = "" )
ggpubr::ggpar(ind.p,title = "")
dev.off()

surv<-as.data.frame(t(adjusted_counts[coef$gene,]))
#PCD_SCORE
coef1<-coef
for (ii in 1:nrow(surv)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*surv[ii,jj]
  }
  surv[ii,20]<-pcd_score
}
colnames(surv)[20]<-"pcd_score"
surv$status<-as.numeric(c(pdata1$`composite_end_point:ch1`,pdata2$`outcome:ch1`,pdata3$outcome,
                          pdata4$`survival status, 0 = censored, 1 = death:ch1`,pdata5$`survival status, 0 = censored, 1 = death:ch1`))
surv$time<-as.numeric(c(pdata1$`survival (months):ch1`,as.numeric(pdata2$time)*12,as.numeric(substr(pdata3$`time to outcome`,1,5))*12,
                        as.numeric(pdata4$`time to death (days):ch1`)/30,as.numeric(pdata5$`time to death (days):ch1`)/30))
surv$age<-as.numeric(c(pdata1$`age (years):ch1`,pdata2$`age:ch1`,pdata3$`age of blood draw`,
                       pdata4$`age:ch1`,pdata5$`age:ch1`))
surv$gender<-c(pdata1$`gender:ch1`,pdata2$`gender:ch1`,pdata3$gender,
               pdata4$`sex (0=female, 1=male):ch1`,pdata5$`sex (0=female, 1=male):ch1`)
surv$gender[which(surv$gender %in% c("1","female"))]<-"Female"
surv$gender[which(surv$gender %in% c("0","male"))]<-"Male"
surv<-surv[complete.cases(surv),]

res.cut <- surv_cutpoint(surv, time = "time", event = "status", variables ="pcd_score")
res.cat <- surv_categorize(res.cut)
surv$group<-res.cat$pcd_score
surv$group[which(surv$group %in% "high")]<-"s_high"
colnames(res.cat)
fit <- survfit(Surv(time,status) ~ pcd_score, data =res.cat)
p<-ggsurvplot(fit, data =res.cat, pval = T,  conf.int = F,risk.table = F, 
              legend.title="pcd_score", title="Overall survival",
              ylab="percentage",xlab = " Time (Months)",
              risk.table.col = "strata", linetype = "strata", 
              surv.median.line = "hv", ggtheme = theme_bw(),  palette = c("#FF6666","#3399FF")) 
pdf(file = "ALLdata_os_surv.pdf",width = 9,height =6.6)
p
dev.off()

res.cut <- surv_cutpoint(surv, time = "time", event = "status", variables =colnames(surv)[1:20])
res.cat <- surv_categorize(res.cut)
res.cat[res.cat == "high"] = "s_high"
univ_formulas <- sapply(colnames(surv)[1:20],function(x) as.formula(paste('Surv(time,status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data =res.cat )})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         HR <-signif(x$coef[2], digits=3);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         res<-c(p.value,HR,HR.confint.lower,HR.confint.upper )
                         names(res)<-c("p.value","HR","HR.confint.lower","HR.confint.upper")
                         return(res)
                       })
res<-data.frame()
for (i in 1:(ncol(surv)-5)) {
  res[i,1]<-names(univ_results)[i]
  res[i,2]<-univ_results[[i]][1]
  res[i,3]<-univ_results[[i]][2]
  res[i,4]<-univ_results[[i]][3]
  res[i,5]<-univ_results[[i]][4]
}
colnames(res)<-c("character","p.value","HR","HR.confint.lower","HR.confint.upper")


#C
library(survival)
library(survminer)
#single_cox
covariates <- colnames(surv)[c(25,23,24)]
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(time,status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = surv)})
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         HR <-signif(x$coef[2], digits=3);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         res<-c(p.value,HR,HR.confint.lower,HR.confint.upper )
                         names(res)<-c("p.value","HR","HR.confint.lower","HR.confint.upper")
                         return(res)
                       })
res<-data.frame()
for (i in 1:3) {
  res[i,1]<-names(univ_results)[i]
  res[i,2]<-univ_results[[i]][1]
  res[i,3]<-univ_results[[i]][2]
  res[i,4]<-univ_results[[i]][3]
  res[i,5]<-univ_results[[i]][4]
}
colnames(res)<-c("character","p.value","HR","HR.confint.lower","HR.confint.upper")
write.table(res,"single_cox.txt",sep = "\t",quote = F,row.names = F)

#multi_cox
s<-Surv(time,status) ~group+age+gender
model<-coxph(s, data =surv )
a<-summary(model,data=surv)[["conf.int"]] 
b<-summary(model,data=surv)[["coefficients"]]
res_multi<-as.data.frame(cbind(a,b))
res_multi$ch<-rownames(res_multi)
res_multi<-res_multi[,c(10,9,1,3,4)]
res_multi[,3:5]<-round(res_multi[,3:5],3)
colnames(res_multi)<-c("character","p.value","HR","HR.confint.lower","HR.confint.upper")
res_multi$p.value<- signif(res_multi$p.value, digits=3)
write.table(res_multi,"multi_cox.txt",sep = "\t",quote = F,row.names = F)


#D
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(utils)
library(ggplot2)
library(corrgram)
library(ggthemes)
library(ggpubr)
#GSE32537
gset=getGEO('GSE32537',AnnotGPL=F)
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
table(pdata$`final diagnosis:ch1`)
pdata<-pdata[which(pdata$`final diagnosis:ch1` %in% c("IPF/UIP")),]
exp<-exp[intersect(rownames(exp),coef$gene),]
surv<-as.data.frame(t(exp[,rownames(pdata)]))
surv$age<-as.numeric(pdata$`age:ch1`)
surv$DLCO<-as.numeric(pdata$`dlco % predicted:ch1`)
surv$FVC<-as.numeric(pdata$`fvc pre-bronchodilator % predicted:ch1`)
coef1<-coef[intersect(rownames(exp),coef$gene),]
for (ii in 1:nrow(surv)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*surv[ii,jj]
  }
  surv[ii,21]<-pcd_score
}
colnames(surv)[21]<-"pcd_score"
p <- ggplot(surv, aes(pcd_score,FVC)) +
  geom_point() +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#9966CC","#FFCC33"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=surv, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="GSE32537-pcd_score_FVC.pdf",width =3,height =3)
p <- ggplot(surv, aes(pcd_score,DLCO)) +
  geom_point() +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#9966CC","#FFCC33"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=surv, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="GSE32537-pcd_score_DLCO.pdf",width =3,height =3)



#E
#GSE47460
gset=getGEO('GSE47460',AnnotGPL=F)
data2 <- data.frame(exprs(gset[[2]]))
expMatrix <- gset[[2]] 
pdata_2<-pData(expMatrix)
gpl2<-gset[["GSE47460-GPL6480_series_matrix.txt.gz"]]@featureData@data
gpl2<-gpl2[rownames(data2),]
exp_2<-aggregate(x=data2[,1:ncol(data2)],by =gpl2$GENE_SYMBOL  %>% list(),FUN = mean)
rownames(exp_2)<-exp_2$Group.1
exp_2<-exp_2[-1,-1]  
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
exp_2<-cbind(exp_2[,pdata_2$geo_accession[c]],exp_2[,pdata_2$geo_accession[d]])
pdata_2<-pdata_2[colnames(exp_2),]
exp<-exp_2[intersect(rownames(exp_2),coef$gene),]
surv<-as.data.frame(t(exp[,rownames(pdata_2)]))
surv$age<-as.numeric(pdata_2$`age:ch1`)
surv$DLCO<-as.numeric(pdata_2$`%predicted dlco:ch1`)
surv$FVC<-as.numeric(pdata_2$`%predicted fvc (pre-bd):ch1`)
surv$FEV1<-as.numeric(pdata_2$`%predicted fev1 (pre-bd):ch1`)
coef1<-coef[intersect(rownames(exp),coef$gene),]
for (ii in 1:nrow(surv)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*surv[ii,jj]
  }
  surv[ii,20]<-pcd_score
}
colnames(surv)[20]<-"pcd_score"
p <- ggplot(surv, aes(pcd_score,FVC)) +
  geom_point() +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#9966CC","#FFCC33"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=surv, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="GSE47460-GPL6480-pcd_score_FVC.pdf",width =3,height =3)
p <- ggplot(surv, aes(pcd_score,DLCO)) +
  geom_point() +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#9966CC","#FFCC33"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=surv, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="GSE47460-GPL6480-pcd_score_DLCO.pdf",width =3,height =3)
p <- ggplot(surv, aes(pcd_score,FEV1)) +
  geom_point() +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#9966CC","#FFCC33"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=surv, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="GSE47460-GPL6480-pcd_score_FEV1.pdf",width =3,height =3)



#F
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
exp_1<-cbind(exp_1[,pdata_1$geo_accession[a]],exp_1[,pdata_1$geo_accession[b]])
pdata_1<-pdata_1[colnames(exp_1),]

exp<-exp_1[intersect(rownames(exp_1),coef$gene),]
surv<-as.data.frame(t(exp[,rownames(pdata_1)]))
surv$age<-as.numeric(pdata_1$`age:ch1`)
surv$DLCO<-as.numeric(pdata_1$`%predicted dlco:ch1`)
surv$FVC<-as.numeric(pdata_1$`%predicted fvc (pre-bd):ch1`)
surv$FEV1<-as.numeric(pdata_1$`%predicted fev1 (pre-bd):ch1`)
coef1<-coef[intersect(rownames(exp),coef$gene),]
for (ii in 1:nrow(surv)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*surv[ii,jj]
  }
  surv[ii,20]<-pcd_score
}
colnames(surv)[20]<-"pcd_score"
p <- ggplot(surv, aes(pcd_score,FVC)) +
  geom_point() +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#9966CC","#FFCC33"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=surv, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="GSE47460-GPL14550-pcd_score_FVC.pdf",width =3,height =3)
p <- ggplot(surv, aes(pcd_score,DLCO)) +
  geom_point() +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#9966CC","#FFCC33"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=surv, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="GSE47460-GPL14550-pcd_score_DLCO.pdf",width =3,height =3)
p <- ggplot(surv, aes(pcd_score,FEV1)) +
  geom_point() +   stat_smooth(method=lm)+          
  theme_bw()+ scale_colour_manual(values=c("#9966CC","#FFCC33"))+
  theme(legend.justification=c(1,0), legend.position=c(1,0),legend.title=element_blank())+
  stat_cor(data=surv, method = "pearson")
gg1<-ggExtra::ggMarginal(p, type = "densigram",
                         xparams=list(fill ="#9966CC"),
                         yparams = list(fill="#FFCC33"))
ggsave(gg1,file="GSE47460-GPL14550-pcd_score_FEV1.pdf",width =3,height =3)




#G
library(GSVA)
library(pheatmap)
library(GSEABase)
library(ggplot2)
library(limma)
geneset <- getGmt(file.path("h.all.v2024.1.Hs.symbols.gmt"))
deg_exp<-as.matrix(data)
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(RColorBrewer)
#GSE175457
gset=getGEO('GSE175457',AnnotGPL=F)
exprset <-as.data.frame(t(read.table("GSE175457_geo_sysbio_rnaseq_data.txt",header =T)))
eg = mapIds(x = org.Hs.eg.db,keys =substr(rownames(exprset),1,15),keytype = "ENSEMBL",column = "SYMBOL")
exp4<-aggregate(x=exprset,by = eg %>% list(),FUN = mean)
rownames(exp4)<-exp4$Group.1
exp4<-exp4[,-1]
pdata4<-gset[["GSE175457_series_matrix.txt.gz"]]@phenoData@data
table(pdata4$characteristics_ch1.3)
gene<-strsplit(pdata4$title," ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])
}
pdata4$title<-gene_symbol
exp4<-exp4[,pdata4$title]
colnames(exp4)<-pdata4$geo_accession
pdata4$status<-pdata4$characteristics_ch1.3
pdata4$status[which(pdata4$status %in% c("case_control_status: control"))]<-"control"
pdata4$status[which(pdata4$status %in% c("case_control_status: case"))]<-"IPF"
exp4<-exp4[,rownames(pdata4)]
ssgsea_par <- ssgseaParam(as.matrix(exp4), geneset)  # all other values are default values
es.max1 <- gsva(ssgsea_par)
rownames(es.max1 )
pathway<-c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_TGF_BETA_SIGNALING",
"HALLMARK_APOPTOSIS","HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY","HALLMARK_WNT_BETA_CATENIN_SIGNALING")

library(ggplot2)  
library(dplyr)  
library(vegan) 
library(linkET)  
library(Hmisc)
library("Rmisc")
library("plyr")
PCD_feature_exp<-as.data.frame(t(data[coef$gene,ipf]))
gsva_hall<-as.data.frame(t(es.max1[pathway,rownames(PCD_feature_exp)]))
colnames(gsva_hall)<-gsub("_"," ",tolower(colnames(gsva_hall)))
immu1<-cbind(PCD_feature_exp,gsva_hall)
cor<-data.frame()
for (mm in 20:24) {
  cor1<-data.frame(0,0,0,0)
  for (i in 1:19) {
    c<-rcorr(immu1[,i],immu1[,mm])
    cor1[i,1]<-colnames(immu1)[mm]
    cor1[i,2]<-colnames(immu1)[i]
    cor1[i,3]<-c$r[2]
    cor1[i,4]<-c$P[2]
  }
  cor<-rbind(cor,cor1)
}
colnames(cor)<-c("spec","env","r","p")  
cor<-cor[which(cor$p<0.01),]
cor$rd<-cut(abs(cor$r), breaks = c(-Inf,0.3, 0.5, Inf), 
            labels = c("< 0.3", "0.3 - 0.5", ">= 0.5"))
cor$pd<-cut(cor$r, breaks = c(-Inf,0,Inf),
            labels = c("< 0", ">0"))
p1<-qcorrplot(correlate(PCD_feature_exp,), type = "upper", diag = FALSE) +  
  geom_square() + 
  geom_couple(aes(xend = .xend + 1.25, 
                  yend = .yend + 0.5, colour = pd, size = rd  ), #
              data = cor,  
              curvature = 0.1 ) +  
  geom_mark(only_mark = TRUE)+  
  geom_diag_label( mapping = aes(y = .y + 0.05), hjust = 0.15) + 
  scale_fill_gradientn(colours = c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7",
                                   "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")) +
  scale_size_manual(values = c(0.5, 1, 2)) + 
  scale_colour_manual(values = color_pal(3)) + 
  guides( size = guide_legend( title = "Mantel's r",    
                               override.aes = list(colour = "grey35"), order = 2),  
          colour = guide_legend( title = "Mantel's p",   
                                 override.aes = list(size = 3),order = 1),       
          fill = guide_colorbar(title = "Pearson's r", order = 3)) +   
  theme( axis.text.y = element_blank() )   
ggsave("PCD_feature-gsva_hall-IPF.pdf", plot= p1, width = 11, height = 8)

PCD_feature_exp<-as.data.frame(t(data[coef$gene,con]))
gsva_hall<-as.data.frame(t(es.max1[pathway,rownames(PCD_feature_exp)]))
colnames(gsva_hall)<-gsub("_"," ",tolower(colnames(gsva_hall)))
immu1<-cbind(PCD_feature_exp,gsva_hall)
cor<-data.frame()
for (mm in 20:24) {
  cor1<-data.frame(0,0,0,0)
  for (i in 1:19) {
    c<-rcorr(immu1[,i],immu1[,mm])
    cor1[i,1]<-colnames(immu1)[mm]
    cor1[i,2]<-colnames(immu1)[i]
    cor1[i,3]<-c$r[2]
    cor1[i,4]<-c$P[2]
  }
  cor<-rbind(cor,cor1)
}
colnames(cor)<-c("spec","env","r","p")  
cor<-cor[which(cor$p<0.01),]
cor$rd<-cut(abs(cor$r), breaks = c(-Inf,0.3, 0.5, Inf), 
            labels = c("< 0.3", "0.3 - 0.5", ">= 0.5"))
cor$pd<-cut(cor$r, breaks = c(-Inf,0,Inf),
            labels = c("< 0", ">0"))
p1<-qcorrplot(correlate(PCD_feature_exp,), type = "lower", diag = FALSE) +  geom_square() +  
  geom_couple(aes(xend = .xend + 1.25,      
                  yend = .yend + 0.5, colour = pd, size = rd  ),
              data = cor,     curvature = 0.1 ) +   
  geom_mark(only_mark = TRUE)+  
  geom_diag_label( mapping = aes(y = .y + 0.05), hjust = 0.15) +  
  scale_fill_gradientn(colours = c("#053061","#2166AC","#4393C3","#92C5DE","#D1E5F0","#F7F7F7",
                                   "#FDDBC7","#F4A582","#D6604D","#B2182B","#67001F")) +  
  scale_size_manual(values = c(0.5, 1, 2)) + 
  scale_colour_manual(values = color_pal(3)) + 
  guides( size = guide_legend( title = "Mantel's r",    
                               override.aes = list(colour = "grey35"), order = 2),  
          colour = guide_legend( title = "Mantel's p",  
                                 override.aes = list(size = 3),order = 1),       
          fill = guide_colorbar(title = "Pearson's r", order = 3)) +   
  theme( axis.text.y = element_blank() )   
ggsave("PCD_feature-gsva_hall-CON.pdf", plot= p1, width = 11, height = 8)




#H
set.seed(7889)
library(pheatmap)
library(ggplot2)
#GSE175457
exp1<-read.table("GSE175457_expr.txt")
pdata1<-read.table("GSE175457_class.txt")
ipf<-pdata1$sample[which(pdata1$outcome %in% "1")]
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
phon<-as.data.frame(t(exp1[coef$gene,pdata1$sample]))
phon$class<-pdata1$outcome
for (ii in 1:nrow(phon)) {
  pcd_score<-0
  for (jj in 1:nrow(coef)) {
    pcd_score<-pcd_score+coef[jj,2]*phon[ii,jj]
  }
  phon[ii,21]<-pcd_score
}
colnames(phon)[21]<-"pcd_score"
data1<-as.data.frame(t(exp1))
data1$pcd_score<-phon$pcd_score
library(Hmisc)
library("Rmisc")
library("plyr")
cor<-data.frame()
for (i in 1:(ncol(data1)-1)) {
  c<-rcorr(data1[,i],data1[,ncol(data1)])
  cor[i,1]<-c$r[2]
  cor[i,2]<-c$P[2]
  cor[i,3]<-colnames(data1)[i]
  cor[i,4]<-colnames(data1)[ncol(data1)]
}
colnames(cor)<-c("correlation","Pvalue","gene","pcd_score")
cor4<-cor[cor$Pvalue<0.01,]
cor4<-cor4[abs(cor4$correlation)>0.7,]
cor4<-na.omit(cor4)
cor4<-cor4[order(cor4$correlation,decreasing = F),]
d-gene-pheatmap.pdf",width = 8,height = 6)

library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
library(enrichplot)
library(ggridges)
library(grid)
library(ggsci)
library(ggplot2)
library(ggthemes)
library(ggpubr)
deg<-cor4$gene[which(cor4$correlation>0)]
eg = mapIds(x = org.Hs.eg.db,keys =deg,keytype = "SYMBOL",column = "ENTREZID")
ego_bp <- enrichGO(OrgDb="org.Hs.eg.db",gene = eg,keyType = "ENTREZID",pvalueCutoff = 0.05,
                   ont = "BP",readable=TRUE)
ego_cc <- enrichGO(OrgDb="org.Hs.eg.db",gene = eg,keyType = "ENTREZID",pvalueCutoff = 0.05,
                   ont = "CC",readable=TRUE)
ego_mf <- enrichGO(OrgDb="org.Hs.eg.db",gene = eg,keyType = "ENTREZID",pvalueCutoff = 0.05,
                   ont = "MF",readable=TRUE)
ego_bp_result<-ego_bp@result[order(ego_bp@result$pvalue),]
ego_cc_result<-ego_cc@result[order(ego_cc@result$pvalue),]
ego_mf_result<-ego_mf@result[order(ego_mf@result$pvalue),]
go_enrich<-rbind(ego_bp_result[1:10,],ego_cc_result[1:10,],ego_mf_result[1:10,])
go_enrich$type<-c(rep("BP",10),rep("CC",10),rep("MF",10))
go_enrich$p<--log10(go_enrich$pvalue)
go_enrich$Description<-factor(go_enrich$Description,levels = go_enrich$Description)
colnames(go_enrich)
p<-ggplot(go_enrich,aes(x=Description,y=p, fill=type)) + #
  geom_bar(stat="identity", width=0.8) +  
  scale_fill_manual(values = c("#CC66FF", "#99CCFF", "#FFCC66") ) + 
  coord_flip() +
  xlab("GO term") + ylab("Gene_Number") + labs(title = "GO Terms Enrich")+theme_bw()
ggsave(p,file="1-PCD-cor-go+.pdf",width = 8,height =10)
write.csv(rbind(ego_bp_result,ego_cc_result,ego_mf_result),"1-PCD-cor-go+.csv")

#I
deg<-cor4$gene[which(cor4$correlation<0)] 
eg = mapIds(x = org.Hs.eg.db,keys =deg,keytype = "SYMBOL",column = "ENTREZID")
ego_bp <- enrichGO(OrgDb="org.Hs.eg.db",gene = eg,keyType = "ENTREZID",pvalueCutoff = 0.05,
                   ont = "BP",readable=TRUE)
ego_cc <- enrichGO(OrgDb="org.Hs.eg.db",gene = eg,keyType = "ENTREZID",pvalueCutoff = 0.05,
                   ont = "CC",readable=TRUE)
ego_mf <- enrichGO(OrgDb="org.Hs.eg.db",gene = eg,keyType = "ENTREZID",pvalueCutoff = 0.05,
                   ont = "MF",readable=TRUE)
ego_bp_result<-ego_bp@result[order(ego_bp@result$pvalue),]
ego_cc_result<-ego_cc@result[order(ego_cc@result$pvalue),]
ego_mf_result<-ego_mf@result[order(ego_mf@result$pvalue),]
go_enrich<-rbind(ego_bp_result[1:10,],ego_cc_result[1:10,],ego_mf_result[1:10,])
go_enrich$type<-c(rep("BP",10),rep("CC",10),rep("MF",10))
go_enrich$p<--log10(go_enrich$pvalue)
go_enrich$Description<-factor(go_enrich$Description,levels = go_enrich$Description)
colnames(go_enrich)
p<-ggplot(go_enrich,aes(x=Description,y=p, fill=type)) +
  geom_bar(stat="identity", width=0.8) + 
  scale_fill_manual(values = c("#CC66FF", "#99CCFF", "#FFCC66") ) +  
  coord_flip() +  
  xlab("GO term") + ylab("Gene_Number") + labs(title = "GO Terms Enrich")+theme_bw()
ggsave(p,file="1-PCD-cor-go-.pdf",width = 8,height =10)
write.csv(rbind(ego_bp_result,ego_cc_result,ego_mf_result),"1-PCD-cor-go-.csv")
