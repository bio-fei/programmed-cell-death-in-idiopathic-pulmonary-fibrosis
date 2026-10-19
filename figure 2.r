#A
library(GEOquery)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)

#GSE175457 training data
gset=getGEO('GSE175457',AnnotGPL=F)
exprset <-as.data.frame(t(read.table("GSE175457_geo_sysbio_rnaseq_data.txt",header =T)))
eg = mapIds(x = org.Hs.eg.db,keys =substr(rownames(exprset),1,15),keytype = "ENSEMBL",column = "SYMBOL")
exp4<-aggregate(x=exprset,by = eg %>% list(),FUN = mean)
rownames(exp4)<-exp4$Group.1
exp4<-exp4[,-1]
pdata4<-gset[["GSE175457_series_matrix.txt.gz"]]@phenoData@data
gene<-strsplit(pdata4$title," ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])
}
pdata4$title<-gene_symbol
exp4<-exp4[,pdata4$title]
colnames(exp4)<-pdata4$geo_accession
pdata4$status<-pdata4$characteristics_ch1.3
pdata4$status[which(pdata4$status %in% c("case_control_status: control"))]<-"0"
pdata4$status[which(pdata4$status %in% c("case_control_status: case"))]<-"1"

PCD_gene<-read.table("PCD_DEG_104.txt")
exp4<-exp4[PCD_gene$V1,rownames(pdata4)]
pdata4<-pdata4[,2:3]
colnames(pdata4)<-c("sample","outcome")
write.table(exp4,"Training_expr.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(pdata4,"Training_class.txt",col.names = T,row.names = F,sep = "\t",quote = F)

#Testing data
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
pdata<-pdata[which(pdata$`final diagnosis:ch1` %in% c("control","IPF/UIP")),]
pdata$status<-pdata$`final diagnosis:ch1` 
pdata$status[which(pdata$status %in% c("IPF/UIP"))]<-"IPF"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene$V1,pdata$geo_accession]
pdata<-data.frame(Sample=pdata$geo_accession,
                  Cohort=rep("GSE32537",nrow(pdata)),
                  outcome=pdata$status)
exp[is.na(exp)]=0
rownames(exp)<-PCD_gene$V1
test_exp<-exp
test_class<-pdata

#GSE47460
gset=getGEO('GSE47460',AnnotGPL=F)
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
exp<-exp[PCD_gene$V1,pdata$geo_accession]
pdata<-data.frame(Sample=pdata$geo_accession,
                  Cohort=rep("GSE47460_GPL14550",nrow(pdata)),
                  outcome=pdata$status)
exp[is.na(exp)]=0
rownames(exp)<-PCD_gene$V1
test_exp<-cbind(test_exp,exp)
test_class<-rbind(test_class,pdata)

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
exp<-cbind(exp_2[,pdata_2$geo_accession[c]],exp_2[,pdata_2$geo_accession[d]])
pdata<-pdata_2
pdata<-pdata[colnames(exp),]
pdata$status<-c(rep("IPF",length(c)),rep("control",length(d)))
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp<-exp[PCD_gene$V1,pdata$geo_accession]
pdata<-data.frame(Sample=pdata$geo_accession,
                  Cohort=rep("GSE47460_GPL6480",nrow(pdata)),
                  outcome=pdata$status)
exp[is.na(exp)]=0
rownames(exp)<-PCD_gene$V1
test_exp<-cbind(test_exp,exp)
test_class<-rbind(test_class,pdata)

write.table(test_exp,"Testing_expr.txt",col.names = T,row.names = T,sep = "\t",quote = F)
write.table(test_class,"Testing_class.txt",col.names = T,row.names = F,sep = "\t",quote = F)


work.path <- "PCDS"
setwd(work.path) 
code.path <- file.path(work.path, "Codes")
data.path <- file.path(work.path, "InputData")
res.path <- file.path(work.path, "Results")
fig.path <- file.path(work.path, "Figures")

library(openxlsx)
library(seqinr)
library(plyr)
library(randomForestSRC)
library(glmnet)
library(plsRglm)
library(gbm)
library(caret)
library(mboost)
library(e1071)
library(BART)
library(MASS)
library(snowfall)
library(xgboost)
library(ComplexHeatmap)
library(RColorBrewer)
library(pROC)

source(file.path(code.path, "ML.R"))
FinalModel <- c("panML", "multiLogistic")[2]

## Training Cohort ---------------------------------------------------------
Train_expr <- read.table(file.path(data.path, "Training_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Train_class <- read.table(file.path(data.path, "Training_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Train_class), colnames(Train_expr))
Train_expr <- Train_expr[,comsam]; Train_class <- Train_class[comsam,,drop = F]

## Validation Cohort -------------------------------------------------------
Test_expr <- read.table(file.path(data.path, "Testing_expr.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
Test_class <- read.table(file.path(data.path, "Testing_class.txt"), header = T, sep = "\t", row.names = 1,check.names = F,stringsAsFactors = F)
comsam <- intersect(rownames(Test_class), colnames(Test_expr))
Test_expr <- Test_expr[,comsam];Test_class <- Test_class[comsam,,drop = F]


comgene <- intersect(rownames(Train_expr),rownames(Test_expr))
Train_expr <- t(Train_expr[comgene,]) 
Test_expr <- t(Test_expr[comgene,]) 

Train_set = scaleData(data = Train_expr, centerFlags = F, scaleFlags = F) 
names(x = split(as.data.frame(Test_expr), f = Test_class$Cohort))
Test_set = scaleData(data = Test_expr, cohort = Test_class$Cohort, centerFlags = F, scaleFlags = F)
methods <- read.xlsx(file.path(code.path, "methods.xlsx"), startRow = 2)
methods <- methods$Model
methods <- gsub("-| ", "", methods)

## Train the model --------------------------------------------------------

classVar = "outcome" 
min.selected.var = 10 
Variable = colnames(Train_set)
preTrain.method =  strsplit(methods, "\\+") 
preTrain.method = lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method = unique(unlist(preTrain.method)) 

preTrain.var <- list() 
set.seed(seed = 777) 
for (method in preTrain.method){
  preTrain.var[[method]] = RunML(method = method, Train_set = Train_set, 
                                 Train_label = Train_class, mode = "Variable",       classVar = classVar) }
preTrain.var[["simple"]] <- colnames(Train_set)

## Model training
model <- list() 
set.seed(seed = 777) 
Train_set_bk = Train_set
for (method in methods){
  cat(match(method, methods), ":", method, "\n")
  method_name = method 
  method <- strsplit(method, "\\+")[[1]] 
  
  if (length(method) == 1) method <- c("simple", method) 
  
  Variable = preTrain.var[[method[1]]] 
  Train_set = Train_set_bk[, Variable]  
  Train_label = Train_class            
  model[[method_name]] <- RunML(method = method[2],        
                                Train_set = Train_set,     
                                Train_label = Train_label, 
                                mode = "Model",            
                                classVar = classVar)       
  if(length(ExtractVar(model[[method_name]])) <= min.selected.var) {
    model[[method_name]] <- NULL
  }
}
Train_set = Train_set_bk; rm(Train_set_bk) 
saveRDS(model, file.path(res.path, "model.rds")) 

if (FinalModel == "multiLogistic"){
  logisticmodel <- lapply(model, function(fit){ 
    tmp <- glm(formula = Train_class[[classVar]] ~ .,
               family = "binomial", 
               data = as.data.frame(Train_set[, ExtractVar(fit)]))
    tmp$subFeature <- ExtractVar(fit)
    return(tmp)
  })
}
saveRDS(logisticmodel, file.path(res.path, "logisticmodel.rds")) 

## Evaluate the model -----------------------------------------------------
model <- readRDS(file.path(res.path, "model.rds"))
methodsValid <- names(model)
model[["glmBoost+Ridge"]]$beta

RS_list <- list()
for (method in methodsValid){
  print(method)
  RS_list[[method]] <- CalPredictScore(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) 
}
RS_mat <- as.data.frame(t(do.call(rbind, RS_list)))
write.table(RS_mat, file.path(res.path, "RS_mat.txt"),sep = "\t", row.names = T, col.names = NA, quote = F) 

Class_list <- list()
for (method in methodsValid){
  Class_list[[method]] <- PredictClass(fit = model[[method]], 
                                       new_data = rbind.data.frame(Train_set,Test_set)) 
}
Class_mat <- as.data.frame(t(do.call(rbind, Class_list)))
write.table(Class_mat, file.path(res.path, "Class_mat.txt"), #
            sep = "\t", row.names = T, col.names = NA, quote = F)

fea_list <- list()
for (method in methodsValid) {
  fea_list[[method]] <- ExtractVar(model[[method]])
}

fea_df <- lapply(model, function(fit){
  data.frame(ExtractVar(fit))
})
fea_df <- do.call(rbind, fea_df)
fea_df$algorithm <- gsub("(.+)\\.(.+$)", "\\1", rownames(fea_df))
colnames(fea_df)[1] <- "features"
write.table(fea_df, file.path(res.path, "fea_df.txt"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

AUC_list <- list()
for (method in methodsValid){
  print(method)
  AUC_list[[method]] <- RunEval(fit = model[[method]],     
                                Test_set = Test_set,      
                                Test_label = Test_class,  
                                Train_set = Train_set,    
                                Train_label = Train_class, 
                                Train_name = "GSE175457",   
                                cohortVar = "Cohort",      
                                classVar = classVar)      
}
AUC_mat <- do.call(rbind, AUC_list)
write.table(AUC_mat, file.path(res.path, "AUC_mat.txt"),
            sep = "\t", row.names = T, col.names = T, quote = F)

# Plot --------------------------------------------------------------------
AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)           
avg_AUC <- sort(avg_AUC, decreasing = T)     
AUC_mat <- AUC_mat[names(avg_AUC), ]      
fea_sel <- fea_list[[rownames(AUC_mat)[1]]] 
avg_AUC <- as.numeric(format(avg_AUC, digits = 7, nsmall = 3)) #

if(ncol(AUC_mat) < 3) { 
  CohortCol <- c("#FF9900","#663399") 
} else { 
  CohortCol <- c(brewer.pal(n = 10, name = "Set3"),brewer.pal(n = 7, name = "Paired")) 
}
names(CohortCol) <- colnames(AUC_mat)

cellwidth = 1; cellheight = 0.5
hm <- SimpleHeatmap(AUC_mat, 
                    avg_AUC, 
                    CohortCol, "#c490a0", 
                    cellwidth = cellwidth, cellheight = cellheight, 
                    cluster_columns = F, cluster_rows = F) 
pdf(file.path(fig.path, "AUC1.pdf"), width = cellwidth * ncol(AUC_mat) + 3, height = cellheight * nrow(AUC_mat) * 0.45)
draw(hm)
invisible(dev.off())

#B
res.path <- file.path(work.path, "Results")
AUC_mat <- read.table(file.path(res.path, "AUC_mat.txt"),sep = "\t", row.names = 1, header = T,check.names = F,stringsAsFactors = F)
avg_AUC <- apply(AUC_mat, 1, mean)           
avg_AUC <- sort(avg_AUC, decreasing = F)     
fea_list<-read.table("fea_df.txt",header = T)
fea_preq<-as.data.frame(table(fea_list$algorithm))
rownames(fea_preq)<-fea_preq$Var1
fea_preq<-fea_preq[names(avg_AUC), ]
fea_preq$Var1<-factor(fea_preq$Var1,levels = fea_preq$Var1)
z<-fea_preq$Freq
p<-ggplot(data = fea_preq, mapping = aes(x = Freq,y=Var1,fill = Freq)) +
  geom_bar(stat = 'identity',colour = 'gray')+
  scale_fill_gradient(low = "#FFFAFA", high = "#FF9900") +
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+
  theme(panel.background=element_rect(fill="white",colour="white"), 
        axis.line=element_line(colour="black"))+
  xlab(" ")+ylab('algorithm')+labs(fill='')+
  geom_text(mapping = aes(label = z))
ggsave(p,file="105_algorithm_fea_freq.pdf",width = 8,height = 15)


model<-readRDS("model.rds")
coef<-as.matrix(model[["glmBoost+Ridge"]]$beta)
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])

data<-read.table("Training_expr.txt")
data_calss<-read.table("Training_class.txt",header = T)
ipf<-data_calss$sample[which(data_calss$outcome %in% "1")]
con<-data_calss$sample[which(data_calss$outcome %in% "0")]

coef1<-coef
coef1<-coef1[order(coef1$coef,decreasing = F),]
coef1$gene<-factor(coef1$gene,levels = coef1$gene)
z<-round(coef1$coef,digits = 2)
p<-ggplot(data = coef1, mapping = aes(x = coef,y=gene,fill = coef)) +
  geom_bar(stat = 'identity',colour = 'gray')+
  scale_fill_gradient2(low = "#6699FF",mid = "white", high = "#CB5746") +
  theme(axis.text.x = element_text(angle = 0,hjust = 1))+
  theme(panel.background=element_rect(fill="white",colour="white"), 
        axis.line=element_line(colour="black"))+
  xlab("Coefficient of variation")+ylab('')+labs(fill='')+
  geom_text(mapping = aes(label = z))
ggsave(p,file="fea_gene19_Coef.pdf",width = 6,height = 6)



#C-D

#GSE175457
model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
data<-read.table("Training_expr.txt")
data_calss<-read.table("FTraining_class.txt",header = T)
ipf<-data_calss$sample[which(data_calss$outcome %in% "1")]
con<-data_calss$sample[which(data_calss$outcome %in% "0")]
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

library(ggDCA)
library(rmda)
library(rms)
phon1<-phon
phon1$class <- phon1$class%>% as.numeric()
colnames(phon1)
pcd_score<-decision_curve(class ~ pcd_score, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
CDH3<-decision_curve(class ~ CDH3, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
MMP11<-decision_curve(class~ MMP11, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
MMP13<-decision_curve(class~ MMP13, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
CDH2<-decision_curve(class ~ CDH2, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
PROM2<-decision_curve(class~ PROM2, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
IL13RA2<-decision_curve(class~ IL13RA2, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
ADRA1A<-decision_curve(class~ADRA1A, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
SFRP2<-decision_curve(class ~SFRP2, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
CD24<-decision_curve(class ~ CD24, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
TWIST1<-decision_curve(class~TWIST1, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
CLDN1<-decision_curve(class ~CLDN1, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
PTHLH<-decision_curve(class ~PTHLH, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
FZD5<-decision_curve(class ~ FZD5, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
EDA2R<-decision_curve(class ~EDA2R, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
CCK<-decision_curve(class ~  CCK, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
NFATC4<-decision_curve(class~NFATC4, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
FAIM<-decision_curve(class ~ FAIM, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
CDCP1<-decision_curve(class ~CDCP1, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")
UVRAG<-decision_curve(class ~UVRAG, data = phon1, family = binomial(link = "logit"),thresholds = seq(0, 1, by = 0.01),confidence.intervals = "none")

pdf(width = 12,height = 8,file = "DAC_score.pdf")
combined_dca <- list("PCD score"=pcd_score,
            "CDH3"=CDH3,"MMP11"=MMP11,"MMP13"=MMP13,"CDH2"=CDH2,"PROM2"=PROM2,"IL13RA2"=IL13RA2,
            "ADRA1A"=ADRA1A,"SFRP2"=SFRP2,"CD24"=CD24,"TWIST1"=TWIST1,"CLDN1"=CLDN1,"PTHLH"=PTHLH,
            "FZD5"=FZD5,"EDA2R"=EDA2R,"CCK"=CCK,"NFATC4"=NFATC4,"FAIM"=FAIM,"CDCP1"=CDCP1,"UVRAG"=UVRAG)
plot_decision_curve(combined_dca,curve.names=names(combined_dca),cost.benefit.axis = FALSE,
                    confidence.intervals = FALSE,standardize = TRUE,legend.position = "topright",
                    col= c(pal_d3("category20")(20)[1:14],pal_locuszoom()(7)[6:7],pal_d3("category20")(20)[15:20]))
dev.off()

library(ResourceSelection)
model_glm <- glm(class~ pcd_score,data = phon1, family = binomial)
p.hoslem <- hoslem.test(model_glm$y, fitted(model_glm), g=10)$p.value
p.hoslem
fit2 <- lrm(class~ pcd_score,data = phon1,x=T,y=T)
cal2 <- calibrate(fit2, method='boot', B=500)
pdf(width = 12,height = 8,file = "class_PCDscore.pdf")
plot(cal2,
     xlim = c(0,1),
     ylim = c(0,1),
     xlab = "Prediced Probability",
     ylab = "Observed Probability",
     cex.lab=1.2, cex.axis=1, cex.main=1.2, cex.sub=0.8,
     legend = FALSE
) 
lines(cal2[,c("predy","calibrated.corrected")], 
      type = 'l', 
      lwd = 3, 
      pch = 16, 
      col = "#2166AC") 
lines(cal2[,c("predy","calibrated.orig")],type="l",pch=16,lwd=3,col="tomato")
abline(0,1,
       lty = 2, 
       lwd = 2, 
       col = "#224444")
legend(0.6,0.2,
       c("Apparent","Bias-corrected","Ideal"), 
       lty = c(2,1,1), 
       lwd = c(2,3,3), 
       col = c("black","#2166AC","tomato"), 
       bty = "n"
)
text(0,0,bquote("Hosmer-Lemeshow "~italic(P)~" = "~.(round(p.hoslem,3))),adj = 0)
dev.off()



#E

model<-readRDS("model.rds")
coef<-data.frame(gene=rownames(as.matrix(model[["glmBoost+Ridge"]]$beta)),
                 coef=as.matrix(model[["glmBoost+Ridge"]]$beta)[,1])
PCD_gene<-coef$gene
model_compare<-data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0)
colnames(model_compare)<-c("dataset","PCDscore_low","PCDscore","PCDscore_high",
                           "PMID37337515_low","PMID37337515","PMID37337515_high",
                           "PMID37334348_low","PMID37334348","PMID37334348_high",
                           "PMID38274787_low","PMID38274787","PMID38274787_high")
PMID37337515<-data.frame(gene=c("IL13RA2","ARG1","IL1R2","MS4A2","PROK2", "ZNF436","PIR"),
                         coef=c("1.3788","-0.1996","-0.3534","0.1868","-0.2646","0.1633","-0.2862"))
rownames(PMID37337515)<-PMID37337515$gene
PMID37334348<-data.frame(gene=c("ASPN","SFRP2","SLCO4A1"),
                         coef=c("0.8534","0.7400","-0.8030"))
rownames(PMID37334348)<-PMID37334348$gene
PMID38274787<-data.frame(gene=c("AGRP","BIRC3","CDA","FAM20C","MT1E", "NELL2","SNCA"," ZNF91"),
                         coef=c("-0.417","-0.296","0.426","0.418","0.086","-0.17","0.731"," -0.198"))
rownames(PMID38274787)<-PMID38274787$gene
model_compare[1,1]<-"GSE175457"
gset=getGEO('GSE175457',AnnotGPL=F)
exprset <-as.data.frame(t(read.table("GSE175457_geo_sysbio_rnaseq_data.txt",header =T)))
eg = mapIds(x = org.Hs.eg.db,keys =substr(rownames(exprset),1,15),keytype = "ENSEMBL",column = "SYMBOL")
exp<-aggregate(x=exprset,by = eg %>% list(),FUN = mean)
rownames(exp)<-exp$Group.1
exp<-exp[,-1]
pdata<-gset[["GSE175457_series_matrix.txt.gz"]]@phenoData@data
gene<-strsplit(pdata$title," ")
gene_symbol=NULL
for(i in 1:length(gene)){
  gene_symbol=c(gene_symbol,gene[[i]][1])
}
length(unique(gene_symbol))
pdata$title<-gene_symbol
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status<-pdata$characteristics_ch1.3
pdata$status[which(pdata$status %in% c("case_control_status: control"))]<-"0"
pdata$status[which(pdata$status %in% c("case_control_status: case"))]<-"1"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[1,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[1,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[1,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[1,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[1,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[1,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[1,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[1,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[1,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[1,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[1,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[1,13]<-round(as.numeric(roc[["ci"]][3]),3)

#GSE32537
model_compare[2,1]<-"GSE32537"
gset=getGEO('GSE32537', AnnotGPL=F)
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
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[2,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[2,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[2,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[2,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[2,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[2,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[2,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[2,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[2,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[2,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[2,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[2,13]<-round(as.numeric(roc[["ci"]][3]),3)


#GSE47460
model_compare[3,1]<-"GSE47460_GPL14550"
gset=getGEO('GSE47460', AnnotGPL=F)
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
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[3,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[3,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[3,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[3,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[3,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[3,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[3,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[3,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[3,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[3,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[3,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[3,13]<-round(as.numeric(roc[["ci"]][3]),3)


model_compare[4,1]<-"GSE47460_GPL6480"
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
exp<-cbind(exp_2[,pdata_2$geo_accession[c]],exp_2[,pdata_2$geo_accession[d]])
pdata<-pdata_2
pdata<-pdata[colnames(exp),]
pdata$status<-c(rep("IPF",length(c)),rep("control",length(d)))
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[4,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[4,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[4,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[4,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[4,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[4,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[4,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[4,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[4,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[4,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[4,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[4,13]<-round(as.numeric(roc[["ci"]][3]),3)


#GSE150910
model_compare[5,1]<-"GSE150910"
gset=getGEO('GSE150910', AnnotGPL=F)
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
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[5,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[5,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[5,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[5,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[5,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[5,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[5,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[5,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[5,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[5,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[5,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[5,13]<-round(as.numeric(roc[["ci"]][3]),3)


 
#GSE134692
model_compare[6,1]<-"GSE134692"
gset=getGEO('GSE134692', AnnotGPL=F)
exprset <-read.table(" GSE134692_tmm_normalized_and_filtered_log2CPM.txt",header = T)
eg = mapIds(x = org.Hs.eg.db,keys =exprset$ensembl,keytype = "ENSEMBL",column = "SYMBOL")
length(unique(eg))
exp<-aggregate(x=exprset[,2:ncol(exprset)],by = eg %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE134692_series_matrix.txt.gz"]]@phenoData@data
pdata<-pdata[which(pdata$characteristics_ch1.2 %in% c("diseasestatus: IPF","diseasestatus: Normal")),]
pdata$status<-pdata$characteristics_ch1.2 
pdata$status[which(pdata$status %in% c("diseasestatus: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diseasestatus: Normal"))]<-"control"
exp<-exp[,pdata$description]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[6,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[6,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[6,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[6,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[6,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[6,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[6,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[6,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[6,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[6,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[6,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[6,13]<-round(as.numeric(roc[["ci"]][3]),3)


#GSE124685
model_compare[7,1]<-"GSE124685"
gset=getGEO('GSE124685', AnnotGPL=F)
exprset <-read.table(" GSE124685_FPKM_Matrix_84_samples.txt",header = T)
exp<-aggregate(x=exprset[,9:ncol(exprset)],by = exprset$gene_id %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-log2(exp[-1,-1]+1)
pdata<-gset[["GSE124685_series_matrix.txt.gz"]]@phenoData@data
pdata$status<-pdata$characteristics_ch1
pdata$status[which(pdata$status %in% c("disease state: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("disease state: healthy donor"))]<-"control"
exp<-exp[,pdata$geo_accession]
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[7,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[7,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[7,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[7,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[7,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[7,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[7,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[7,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[7,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[7,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[7,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[7,13]<-round(as.numeric(roc[["ci"]][3]),3)


 
#GSE110147
model_compare[8,1]<-"GSE110147"
gset=getGEO('GSE110147', AnnotGPL=F)
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
pdata<-pdata[which(pdata$characteristics_ch1 %in% c("disease state: Idiopathic pulmonary fibrosis",
                                                    "disease state: Normal control")),]
pdata$status<-pdata$characteristics_ch1
pdata$status[which(pdata$status %in% c("disease state: Idiopathic pulmonary fibrosis"))]<-"IPF"
pdata$status[which(pdata$status %in% c("disease state: Normal control"))]<-"control"
exp<-exp[,pdata$geo_accession]
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[8,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[8,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[8,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[8,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[8,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[8,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[8,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[8,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[8,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[8,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[8,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[8,13]<-round(as.numeric(roc[["ci"]][3]),3)


 
#GSE53845
model_compare[9,1]<-"GSE53845"
gset=getGEO('GSE53845', AnnotGPL=F)
exprset <- data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE53845_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE53845_series_matrix.txt.gz"]]@phenoData@data
pdata$status<-pdata$characteristics_ch1.2
pdata$status[which(pdata$status %in% c("diagnosis: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diagnosis: control"))]<-"control"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[9,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[9,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[9,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[9,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[9,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[9,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[9,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[9,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[9,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[9,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[9,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[9,13]<-round(as.numeric(roc[["ci"]][3]),3)


 
#GSE92592
model_compare[10,1]<-"GSE92592"
gset=getGEO('GSE92592', AnnotGPL=F)
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
pdata$status<-pdata$`disease state:ch1`
pdata$status[which(pdata$status %in% c("Idiopathic Pulmonary Fibrosis"))]<-"IPF"
pdata$status[which(pdata$status %in% c("Control"))]<-"control"
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[10,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[10,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[10,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[10,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[10,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[10,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[10,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[10,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[10,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[10,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[10,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[10,13]<-round(as.numeric(roc[["ci"]][3]),3)


 
#GSE231693
model_compare[11,1]<-"GSE231693"
gset=getGEO('GSE231693', AnnotGPL=F)
exprset <- read.table(" GSE231693_count.tsv",header =T)
eg = mapIds(x = org.Hs.eg.db,keys =rownames(exprset),keytype = "ENSEMBL",column = "SYMBOL")
length(unique(eg))
exp<-aggregate(x=exprset,by = eg %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-log2(exp[,-1]+1)
pdata<-gset[["GSE231693_series_matrix.txt.gz"]]@phenoData@data
pdata<-pdata[which(pdata$characteristics_ch1.1 %in% c("disease state: idiopathic pulmonary fibrosis","disease state: normal")),]
pdata$status<-pdata$characteristics_ch1.1
pdata$status[which(pdata$status %in% c("disease state: idiopathic pulmonary fibrosis"))]<-"IPF"
pdata$status[which(pdata$status %in% c("disease state: normal"))]<-"control"
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[11,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[11,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[11,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[11,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[11,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[11,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[11,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[11,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[11,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[11,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[11,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[11,13]<-round(as.numeric(roc[["ci"]][3]),3)


 
#GSE213001
model_compare[12,1]<-"GSE213001"
gset=getGEO('GSE213001', AnnotGPL=F)
exprset <- as.data.frame(fread(" GSE213001_Entrez-IDs-Lung-IPF-GRCh38-p12-raw_counts.csv",header =T))
eg = mapIds(x = org.Hs.eg.db,keys =exprset$V1,keytype = "ENSEMBL",column = "SYMBOL")
length(unique(eg))
exp<-aggregate(x=exprset[,-1],by = eg %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-log2(exp[,-1]+1)
pdata<-gset[["GSE213001_series_matrix.txt.gz"]]@phenoData@data
pdata<-pdata[which(pdata$characteristics_ch1.10 %in% c("diagnosis: IPF","diagnosis: NDC")),]
pdata$status<-pdata$characteristics_ch1.10
pdata$status[which(pdata$status %in% c("diagnosis: IPF"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diagnosis: NDC"))]<-"control"
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[12,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[12,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[12,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[12,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[12,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[12,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[12,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[12,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[12,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[12,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[12,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[12,13]<-round(as.numeric(roc[["ci"]][3]),3)



 
#GSE173355
model_compare[13,1]<-"GSE173355"
gset=getGEO('GSE173355', AnnotGPL=F)
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
pdata$status<-pdata$characteristics_ch1.4
pdata$status[which(pdata$status %in% c("diagnosis: ipf"))]<-"IPF"
pdata$status[which(pdata$status %in% c("diagnosis: normal"))]<-"control"
pdata$title<-substr(pdata$title,16,100)
pdata$title<-gsub("-","\\.",pdata$title)
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[13,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[13,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[13,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[13,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[13,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[13,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[13,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[13,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[13,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[13,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[13,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[13,13]<-round(as.numeric(roc[["ci"]][3]),3)



 
#GSE99621
model_compare[14,1]<-"GSE99621"
gset=getGEO('GSE99621', AnnotGPL=F)
exprset <-read.table(" GSE99621_datamatrix.txt",header =T)
rownames(exprset)<-exprset$Symbol
exp<-log2(exprset[,-1]+1)
pdata<-gset[["GSE99621_series_matrix.txt.gz"]]@phenoData@data
pdata$status<-pdata$characteristics_ch1.1
pdata$status[which(pdata$status %in% c("disease_status: idiopathic pulmonary fibrosis"))]<-"IPF"
pdata$status[which(pdata$status %in% c("disease_status: healthy control"))]<-"control"
exp<-exp[,pdata$description.1]
colnames(exp)<-pdata$geo_accession
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[14,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[14,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[14,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[14,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[14,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[14,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[14,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[14,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[14,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[14,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[14,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[14,13]<-round(as.numeric(roc[["ci"]][3]),3)


 
#GSE199949
model_compare[15,1]<-"GSE199949"
gset=getGEO('GSE199949', AnnotGPL=F)
exprset <-as.data.frame(fread(" GSE199949_normalized.txt",header =T))
rownames(exprset)<-exprset$SymboP
exprset <-exprset[,10:51]
exp<-as.data.frame(log2(exprset+1))
pdata<-gset[["GSE199949_series_matrix.txt.gz"]]@phenoData@data
exp<-exp[,pdata$title]
colnames(exp)<-pdata$geo_accession
pdata$status<-pdata$`diagnosis:ch1`
pdata$status[which(pdata$status %in% c("Donor"))]<-"control"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[15,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[15,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[15,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[15,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[15,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[15,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[15,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[15,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[15,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[15,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[15,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[15,13]<-round(as.numeric(roc[["ci"]][3]),3)



#GSE70866
model_compare[16,1]<-"GSE70866"
gset=getGEO('GSE70866', AnnotGPL=F)
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE70866-GPL14550_series_matrix.txt.gz"]]@featureData@data
anno<-anno[rownames(exprset),]
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE70866-GPL14550_series_matrix.txt.gz"]]@phenoData@data
exp<-exp[,pdata$geo_accession]
pdata$status<-pdata$`diagnosis:ch1`
pdata$status[which(pdata$status %in% c("healthy donor"))]<-"control"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[16,2]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[16,3]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[16,4]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[16,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[16,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[16,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[16,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[16,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[16,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[16,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[16,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[16,13]<-round(as.numeric(roc[["ci"]][3]),3)



#GSE33566
model_compare[17,1]<-"GSE33566"
gset=getGEO('GSE33566', AnnotGPL=F) 
exprset <-data.frame(exprs(gset[[1]])) 
anno<-gset[["GSE33566_series_matrix.txt.gz"]]@featureData@data
exp=aggregate(x=exprset,by = anno$GENE_SYMBOL %>% list(),FUN = mean) 
rownames(exp)<-exp$Group.1
exp<-exp[-1,-1]
pdata<-gset[["GSE33566_series_matrix.txt.gz"]]@phenoData@data
exp<-exp[,pdata$geo_accession]
pdata$status<-pdata$source_name_ch1
pdata$status[which(pdata$status %in% c("healthy control (N)_PBB"))]<-"control"
pdata$status[which(pdata$status %in% c("idiopathic pulmonary fibrosis (IPF)_PBB"))]<-"IPF"
pdata$status[which(pdata$status %in% c("IPF"))]<-"1"
pdata$status[which(pdata$status %in% c("control"))]<-"0"
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PCD_gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-coef[intersect(rownames(exp), PCD_gene),]
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}

i<-unique(Test_class$Cohort)[16]
test_auc<-RS_mat[rownames(Test_class)[which(Test_class$Cohort %in% i)],c(1,6)]
test_auc$outcome<-Test_class$outcome[which(Test_class$Cohort %in% i)]
roc <- roc(test_auc$outcome,test_auc$`glmBoost+Ridge` %>% as.numeric(),ci=T,smooth = F)
model_compare[17,2]<-1
model_compare[17,3]<-1
model_compare[17,4]<-1



#PMID37337515
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID37337515$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37337515[intersect(rownames(exp), PMID37337515$gene),]     
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--2.7571
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[17,5]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[17,6]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[17,7]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID37334348
exp_test<-as.data.frame(t(exp[intersect(rownames(exp),PMID37334348$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID37334348[intersect(rownames(exp), PMID37334348$gene),]  
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<--5.6339
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[17,8]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[17,9]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[17,10]<-round(as.numeric(roc[["ci"]][3]),3)

#PMID38274787
exp_test<-as.data.frame(t(exp[intersect(rownames(exp), PMID38274787$gene),rownames(pdata)]))
exp_test$cla<-pdata$status
coef1<-PMID38274787[intersect(rownames(exp), PMID38274787$gene),] 
coef1$coef<-as.numeric(coef1$coef)
for (ii in 1:nrow(exp_test)) {
  pcd_score<-0
  for (jj in 1:nrow(coef1)) {
    pcd_score<-pcd_score+coef1[jj,2]*exp_test[ii,jj]
  }
  exp_test[ii,nrow(coef1)+2]<-pcd_score
}
roc <- roc(exp_test$cla,exp_test[,nrow(coef1)+2] %>% as.numeric(),ci=T,smooth = F)
model_compare[17,11]<-round(as.numeric(roc[["ci"]][1]),3)
model_compare[17,12]<-round(as.numeric(roc[["ci"]][2]),3)
model_compare[17,13]<-round(as.numeric(roc[["ci"]][3]),3)

write.table(model_compare, "model_compare.txt", 
            sep = "\t", row.names = F, col.names = T, quote = F)

rownames(model_compare)<-model_compare$dataset
model_compare<-model_compare[,-1]
model_compare<-as.data.frame(t(model_compare))
model_compare<-model_compare[c(2,5,8,11),]
library(ComplexHeatmap)
library(RColorBrewer)
CohortCol <- c(brewer.pal(n = 10, name = "Set3"),brewer.pal(n = 7, name = "Paired")) 
names(CohortCol) <- colnames(model_compare)
col_ha = columnAnnotation("Cohort" = colnames(model_compare),
                          col = list("Cohort" = CohortCol),
                          show_annotation_name = F)
hm <-Heatmap(as.matrix(model_compare), top_annotation = col_ha,
             col = c( "#FFFFFF", "#CB5746"),cluster_rows = F,cluster_columns = F, 
             rect_gp = gpar(col = "black", lwd = 1), # ?߿�
             show_column_names = T, 
             show_row_names = TRUE,
             row_names_side = "left",
             cell_fun = function(j, i, x, y, w, h, col) { # add text to each grid
               grid.text(label = format(model_compare[i, j], digits = 3, nsmall = 3),
                         x, y, gp = gpar(fontsize = 10))
             })
pdf("model_compare_AUC1.pdf", width = 15, height =5)
draw(hm)
dev.off() 