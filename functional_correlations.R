############Script for correlating functional data from CCLE cell lines#########
## This script contains the code used for comparing cell lines profiles######### 
## using RNA-seq data and methylation data                             #########
################################################################################

####RNA-seq data#####
#Load data
CCLE_expresion<-read.csv("CCLE_expression.csv",header = TRUE) #Load RNA-seq data from DEPMAP
#Link for downloading:https://depmap.org/portal/download/
CCL_names<-read.csv("sample_info (1).csv",header=TRUE) #Load info data with CCLE IDs

#Change Depmap IDs to CCLE IDs
for (i in 1:nrow(CCLE_expresion)){
  CCLE_expresion$X[i]<-CCL_names$cell_line_name[which(CCLE_expresion$X[i]==CCL_names$DepMap_ID)]
}

##Calculate the correlation between OVKATE and ALL cell lines
all<-CCLE_expresion[,-1] #Prepare data for correlations
all<-t(all) #Each column have the data of a cell line
all<-as.data.frame(all)
colnames(all)<-CCLE_expresion$X #Name the colums with CCLE cell lines ID

#Select genes data with abs(logFC)>1 in OVKATE
selected_genes<-c()
for (i in 1:nrow(all)){
  if (abs(all$OVKATE[i])>1){selected_genes<-rbind(selected_genes,all[i,])}
  print(paste0("Linea ",i))
}
saveRDS(selected_genes,"genes_logfc_higher1_ovkate.rds")

#Perform the correlations of gene expression profiles of OVKATE vs all for 
#only the selected genes

cor_all<-c()
for (j in 1:ncol(selected_genes)){
  cor <- cor.test(selected_genes$OVKATE ,selected_genes[,j], method = "pearson")
  cor <- rbind(cellLine_1="OVKATE",cellLine_2=colnames(selected_genes)[j],cor=cor$estimate,p.val=cor$p.value) #sacar rho, p-value, nombre de cell line en una matrix
  cor_all <- cbind(cor_all, cor)
  print(paste0("Comparing n#  ", i))
}
cor_all<-as.data.frame(cor_all)
cor_coef<-as.data.frame(cor_all[3,])
cor_coef<-t(cor_coef)
cor_coef<-as.data.frame(cor_coef)

##Plot correlations of OVKATE vs ALL CELL LINES
ggplot(data=cor_coef,aes(x=as.numeric(cor)))+
  geom_density(fill='firebrick4')+
  xlab("Pearson correlation coefficient of OVKATE vs All cell lines")+
  ylab("Frequency")+
  ggtitle("Distribution of the different correlation coeficients (Gene-expression)")+
  theme_minimal()

# Select correlation results for PANC0203, KURAMOCHI and OV-90
cor_all[which(cor_all[2,]=="OVISE")]
cor_all[which(cor_all[2,]=="KURAMOCHI")]
cor_all[which(cor_all[2,]=="OV-90")]

##Obtention of the  empirical p value
values<-cor_coef[which(cor_coef$cor>0.738),]
empirical_p<-length(values)/nrow(cor_coef)






#### Methylation data###
CCLE_methylation<-read.table("CCLE_RRBS_TSS1kb_20181022.txt",header=T,sep="\t")#Load data from DepMap
methylation_data<-CCLE_methylation[,-c(1:3)] #Prepare data for correlation

#Perform Pearson Correlation except for cell lines whose DNA methylation data 
#was NA
cor_all<-c()
for (j in 1:nrow(methylation_data)){
  if (j==793|j==814|j==822|j==831|j==834|j==835|j==836|j==837|j==838|j==839){}
  else{
  cor <- cor.test(as.numeric(methylation_data$OVKATE_OVARY),as.numeric(methylation_data[,j]), method = "pearson")
  cor <- rbind(cellLine_1="OVKATE",cellLine_2=colnames(methylation_data)[j],cor=cor$estimate) #sacar rho, p-value, nombre de cell line en una matrix
  cor_all <- cbind(cor_all, cor)
  print(paste0("Comparing n#  ", j))}
}
cor_all<-as.data.frame(cor_all)
cor_coef<-as.data.frame(cor_all[3,])
colnames(cor_all)<-c(1:833)
cor_coef<-t(cor_coef)
cor_coef<-as.data.frame(cor_coef)

##Plot correlations of OVKATE vs ALL CELL LINES
ggplot(data=cor_coef,aes(x=as.numeric(cor)))+
  geom_density(fill='firebrick4')+
  xlab("Pearson correlation coefficient of OVKATE vs All cell lines")+
  ylab("Frequency")+
  ggtitle("Distribution of the different correlation coeficients (Gene-expression)")+
  theme_minimal()

# Select correlation results for PANC0203 and OV-90 (There are not data for KUR)
cor_all[which(cor_all[2,]=="PANC0203_PANCREAS")]
cor_all[which(cor_all[2,]=="OV90_OVARY")]

##Obtention of the  empirical p value
values<-cor_coef[which(cor_coef$cor>0.79),]
empirical_p<-length(values)/nrow(cor_coef)

