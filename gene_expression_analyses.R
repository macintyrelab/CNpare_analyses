#####Script for correlating functional data from CCLE and GDSC cell lines#######
## This script contains the code used for comparing cell line profiles ######### 
##  using RNA-seq data                                                 #########
################################################################################

#Load the datasets with Sanger and Broad
broad_data<-read.csv("rnaseq_broad_20210317.csv")
sanger_data<-read.csv("rnaseq_sanger_20210316.csv")
model_id<-read.csv("model_list_20210719.csv") #Data with cell line information
list_broad<-as.data.frame(unique(broad_data$model_id)) #List the cell lines IDs in Broad
list_sanger<-as.data.frame(unique(sanger_data$model_id))#List the cell lines IDs in Sanger

#Select the cell lines which are in both datasets
pairs<-c()
for (i in 1:nrow(list_broad)){
  for (j in 1:nrow(list_sanger)){
    if(list_broad[i,1]==list_sanger[j,1]){
      pairs<-rbind(pairs,list_broad[i,1])
    }
  }
}
#73 pairs of cell line models

#Generate the matrix with comparable data
broad_pairs<-c()
sanger_pairs<-c()
for (i in 1:nrow(pairs)){
  d<-broad_data[which(pairs[i,1]==broad_data$model_id),]
  broad_pairs<-cbind(broad_pairs,d$fpkm)
  
  e<-sanger_data[which(pairs[i,1]==sanger_data$model_id),]
  sanger_pairs<-cbind(sanger_pairs,e$fpkm)
}

#Put names to the columns
broad_pairs<-as.data.frame(broad_pairs)
colnames(broad_pairs)<-pairs
sanger_pairs<-as.data.frame(sanger_pairs)
colnames(sanger_pairs)<-pairs

#Perform correlations
correlations<-getSimilarities(broad_pairs,sanger_pairs,method="pearson")

#Get top hit
colnames(correlations)[3]<-"r"#Needed by getTopHit function to work
list.tophits<-getTopHit(samples=pairs[,1], measure=correlations, method="pearson")

##Determine the percentage of cell lines that are equal
#Top-1 in all measures --> equal profile with equal ploidy
list.tophits$equal<-list.tophits[,1]==list.tophits[,2]
top1<-list.tophits[list.tophits$equal,] 
print("Number of samples in the Top-1 for all comparison measures") 
print(paste0(nrow(top1),"/73 samples"))

