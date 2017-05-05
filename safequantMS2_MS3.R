# lumos and QEP safequant all MS2 quant
# compares two safequeant output files for both MS2 and MS3 TMT quantitation 
setwd("....")

ms2_sq <- read.table("SQ_Results_PROTEIN_MS2.tsv",sep="\t",header=TRUE)
ms3_sq <- read.table("SQ_Results_PROTEIN_MS3.tsv",sep="\t",header=TRUE)

#merge table

sq_merged <- merge(ms2_sq,ms3_sq,all=TRUE,by="proteinDescription")
write.csv(sq_merged,file="sq.merged.csv")


sq_merged.sort <- sq_merged[ , order(names(sq_merged))]
sq_merged.q <- subset(sq_merged.sort,qValue_Condition2.x <= 0.05 | qValue_Condition2.y <= 0.05)
write.csv(sq_merged.q,file="sq_merged.q.csv")
sq_merged.2 <- subset(sq_merged.q,nbPeptides.x >= 2 | nbPeptides.y >= 2)

#calculate % differnt for each column and intensity means

pcchange=function(x,y) (abs(abs(x-y)/(x)))*100

sq_merged.2$log2_diff <- pcchange(sq_merged.2$log2ratio_Condition2.x,sq_merged.2$log2ratio_Condition2.y)
sq_merged.2$pvalue_diff <- pcchange(sq_merged.2$pValue_Condition2.x,sq_merged.2$pValue_Condition2.y)
sq_merged.2$qvalue_diff <- pcchange(sq_merged.2$qValue_Condition2.x,sq_merged.2$qValue_Condition2.y)
sq_merged.2$mean_ms2 <- rowMeans(sq_merged.2[,c(20,22,24,26,28,30,32,34,36,38)])
sq_merged.2$mean_ms3 <- rowMeans(sq_merged.2[,c(21,24,25,27,29,31,33,35,37,9)])

write.csv(sq_merged.2,file="sq.merged.2.csv")

# count number of MS2 proteins without a corresponding MS3


ms2_only_percent <- ((sum(is.na(sq_merged.2$mean_ms3)))/nrow(sq_merged.2))*100
ms3_only_percent <- ((sum(is.na(sq_merged.2$mean_ms2)))/nrow(sq_merged.2))*100

#sort by mean MS2 signal. Theory is that these are the most likly to have interference 

sq_merged.2.sort.ms2 <- sq_merged.2[order(sq_merged.2$mean_ms2),]

# Bottom  and top 10 % of MS2 data

sq_merged.2.sort.ms2.no.na <-sq_merged.2.sort.ms2[rowSums(is.na(sq_merged.2.sort.ms2[,c(20,22,24,26,28,30,32,34,36,38)]))==0,]
nn<- round(nrow(sq_merged.2.sort.ms2.no.na)*0.1)
ms2_bottom10 <-head(sq_merged.2.sort.ms2.no.na, n = nn)
write.csv(ms2_bottom10,file="ms2_bottom10.csv")
ms2_top10 <-tail(sq_merged.2.sort.ms2.no.na, n = nn)
write.csv(ms2_top10,file="ms2_top10.csv")

# Number of MS2 only in top and bottom 10%

ms2_only_percent_top10 <- ((sum(is.na(ms2_top10$mean_ms3)))/nrow(ms2_top10))*100
ms2_only_percent_bottom10 <- ((sum(is.na(ms2_bottom10$mean_ms3)))/nrow(ms2_bottom10))*100


  
#sort by % log_diff

sq_merged.2.sort <- sq_merged.2[order(-sq_merged.2$log2_diff),]
write.csv(sq_merged.2.sort,file="sq_merged.2.sort.csv")




sq_merged.2.sort <- sq_merged.2.sort[!duplicated(sq_merged.2.sort$proteinDescription), ]


 

sq_merged.2.sort.10.1.5fold <- subset(sq_merged.2.sort ,log2ratio_Condition2.x <= -.6 | log2ratio_Condition2.x >= .6)
write.csv(sq_merged.2.sort.10.1.5fold,file="sq_merged.2pep.sort.10percentdifflog2.csv")

oldnames <-colnames(sq_merged.2.sort)
newnames <-gsub(".x",".MS2",oldnames)
newnames <-gsub(".y",".MS3",newnames)
newnames2 <-gsub("X","TMT",newnames)
library(data.table)
setnames(sq_merged.2.sort, oldnames, newnames2)
write.csv(sq_merged.2.sort,file="sq_merged.2.sort.csv")

library(data.table)
sq_merged.2.sort.long <- melt(sq_merged.2.sort, id.vars="proteinDescription",
                                 measure.vars= c("TMT1.MS2","TMT2.MS2",
                                                 "TMT3.MS2","TMT4.MS2",
                                                 "TMT5.MS2","TMT6.MS2",
                                                 "TMT7.MS2","TMT8.MS2",
                                                 "TMT9.MS2","TMT10.MS2",
                                                 "TMT1.MS3","TMT2.MS3",
                                                 "TMT3.MS3","TMT4.MS3",
                                                 "TMT5.MS3","TMT6.MS3",
                                                 "TMT7.MS3","TMT8.MS3",
                                                 "TMT9.MS3","TMT10.MS3"))




group.colors <- c(TMT1.MS2 = "#C0392B",
                  TMT2.MS2 = "#C0392B",
                  TMT3.MS2 = "#C0392B",
                  TMT4.MS2 = "#C0392B",
                  TMT5.MS2 = "#C0392B",
                  TMT6.MS2 = "#E74C3C",
                  TMT7.MS2 = "#E74C3C",
                  TMT8.MS2 = "#E74C3C",
                  TMT9.MS2 = "#E74C3C",
                  TMT10.MS2 = "#E74C3C",
                  TMT1.MS3 = "#2980B9",
                  TMT2.MS3 = "#2980B9",
                  TMT3.MS3 = "#2980B9",
                  TMT4.MS3 = "#2980B9",
                  TMT5.MS3 = "#2980B9",
                  TMT6.MS3 = "#3498DB",
                  TMT7.MS3 = "#3498DB",
                  TMT8.MS3 = "#3498DB",
                  TMT9.MS3 = "#3498DB",
                  TMT10.MS3 = "#3498DB")
                  

group.colors2 <- c(TMT1.MS2 = "#C0392B",
                  TMT2.MS2 = "#C0392B",
                  TMT3.MS2 = "#C0392B",
                  TMT4.MS2 = "#C0392B",
                  TMT5.MS2 = "#C0392B",
                  TMT6.MS2 = "#E74C3C",
                  TMT7.MS2 = "#E74C3C",
                  TMT8.MS2 = "#E74C3C",
                  TMT9.MS2 = "#E74C3C",
                  TMT10.MS2 = "#E74C3C")
        


library(ggplot2)

# log2 plots one pdf


pdf("logplots.pdf")
for (i in 1:nrow(sq_merged.2.sort)) { 
  protein.name <-sq_merged.2.sort[i,15]
  protein.desc <-sq_merged.2.sort[i,13]
  if(is.na(protein.name) == TRUE){
    protein.name=sq_merged.2.sort[i,14]
  }
  test <- subset(sq_merged.2.sort.long,proteinDescription==protein.desc)
  test<- test[complete.cases(test),]
  test$value <-log2(test$value)
  difflog <- paste("%DiffLog2=",round(sq_merged.2.sort[i,32],digits=3))
  log2_MS2 <-paste("log2_MS2=",round(sq_merged.2.sort[i,5],digits=3))
  log2_MS3 <-paste("log2_MS3=",round(sq_merged.2.sort[i,6],digits=3))
  qval_ms2 <-paste("qval_MS2=",round(sq_merged.2.sort[i,18],digits=5))
  qval_ms3 <- paste("qval_MS3=",round(sq_merged.2.sort[i,19],digits=5))
  peps_ms2 <-paste("PepsMS2=",round(sq_merged.2.sort[i,11],digits=3))
  peps_ms3 <- paste("PepsMS3=",round(sq_merged.2.sort[i,12],digits=3))
  cv_ms2 <-  paste("CV_MS2=",round(mean(sq_merged.2.sort[i,1],sq_merged.2.sort[i,3],na.rm = TRUE),digits=3))
  cv_ms3 <-paste("CV_MS3=",round(mean(sq_merged.2.sort[i,2],sq_merged.2.sort[i,4],na.rm = TRUE),digits=3))
  label<- paste(strtrim(protein.name,40), "\n",difflog, "\n",log2_MS2, "\n",log2_MS3, "\n",qval_ms2,"\n",qval_ms3,"\n",peps_ms2,"\n",peps_ms3,"\n",cv_ms2,"\n",cv_ms3)
  if(nrow(test)>10){
    print(ggplot(test, aes(x = variable, y = value)) +
            geom_bar(stat="identity",fill=group.colors)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            geom_text(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=label)))
  } else {
    print(ggplot(test, aes(x = variable, y = value)) +
            geom_bar(stat="identity",fill=group.colors2)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            geom_text(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=label)))
    
  }
}
dev.off()

# regular intensity plots plots one pdf


pdf("plots.pdf")
for (i in 1:nrow(sq_merged.2.sort)) { 
  protein.name <-sq_merged.2.sort[i,15]
  protein.desc <-sq_merged.2.sort[i,13]
  if(is.na(protein.name) == TRUE){
    protein.name=sq_merged.2.sort[i,14]
  }
  test <- subset(sq_merged.2.sort.long,proteinDescription==protein.desc)
  test<- test[complete.cases(test),]
  difflog <- paste("%DiffLog2=",round(sq_merged.2.sort[i,32],digits=3))
  log2_MS2 <-paste("log2_MS2=",round(sq_merged.2.sort[i,5],digits=3))
  log2_MS3 <-paste("log2_MS3=",round(sq_merged.2.sort[i,6],digits=3))
  qval_ms2 <-paste("qval_MS2=",round(sq_merged.2.sort[i,18],digits=5))
  qval_ms3 <- paste("qval_MS3=",round(sq_merged.2.sort[i,19],digits=5))
  peps_ms2 <-paste("PepsMS2=",round(sq_merged.2.sort[i,11],digits=3))
  peps_ms3 <- paste("PepsMS3=",round(sq_merged.2.sort[i,12],digits=3))
  cv_ms2 <-  paste("CV_MS2=",round(mean(sq_merged.2.sort[i,1],sq_merged.2.sort[i,3],na.rm = TRUE),digits=3))
  cv_ms3 <-paste("CV_MS3=",round(mean(sq_merged.2.sort[i,2],sq_merged.2.sort[i,4],na.rm = TRUE),digits=3))
  label<- paste(strtrim(protein.name,40), "\n",difflog, "\n",log2_MS2, "\n",log2_MS3, "\n",qval_ms2,"\n",qval_ms3,"\n",peps_ms2,"\n",peps_ms3,"\n",cv_ms2,"\n",cv_ms3)
  if(nrow(test)>10){
    print(ggplot(test, aes(x = variable, y = value)) +
            geom_bar(stat="identity",fill=group.colors)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            geom_text(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=label)))
  } else {
    print(ggplot(test, aes(x = variable, y = value)) +
            geom_bar(stat="identity",fill=group.colors2)+
            theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            geom_text(aes(x=Inf,y=Inf,hjust=1,vjust=1,label=label)))
    
  }
}
dev.off()
