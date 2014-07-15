
# install the developer version of bioconductor to load the most recent version of cleaver
source("http://bioconductor.org/biocLite.R")
biocLite()
library(BiocInstaller) 
useDevel()
biocLite("cleaver")

#install seqinr
install.packages("seqinr")

library("cleaver")
library("seqinr")

#MacWD
setwd("/Volumes/Lab$/service/on_campus/DePeters/tic analysis/")


#PC WD
#setwd("P:/service/on_campus/DePeters/")

infile <- "Samples View Report for ChineseGoat_Final_050814_tic.xls"
infasta <- "Subset.Database.tic.fasta"

strs<-readLines(infile)
raw.data <- read.table(file=infile, sep="\t",header=TRUE,nrows=length(strs)-1)


sequences <- read.fasta(file=infasta, seqtype="AA", as.string=TRUE, seqonly=FALSE)


sequences <- unlist(sequences)
peptide <- cleave(sequences, custom="[KR]", missedCleavages=0)

Accession.Number<- names(sequences)

# function to calculate number of trypic peptides between 7-30
IBAQ_peps <-function(x) {
  num_aa <- nchar(x)
  ibaq_peps<- x[(num_aa>=7) & (num_aa<=30) ]
  num_ibaq_peps <- length(ibaq_peps)
  return(num_ibaq_peps)
}

ibaqs <-lapply(peptide,IBAQ_peps)

ibaqs2 <-unlist(ibaqs)
ibaqs2<-unname(ibaqs2)
ibaqs_all <-data.frame(Accession.Number,ibaqs2)

# write file with IBAQ peptide numbers per protein
write.csv(ibaqs_all, file="ibaqs_final_tic.csv")

#merge df's
 merged<- merge(ibaqs_all,raw.data,by="Accession.Number")
merged$Protein.Grouping.Ambiguity <-NULL
merged$Quantitative.Variance <-NULL
merged$Taxonomy <-NULL
ibaq_nums<-merged[,1:4]

#divide total peak are by IBAQ numbers
for (i in 1:nrow(merged)) {
  for (j in 5:22)
  ibaq_nums[i,j] <- merged[i,j]/merged[i,2]
}

goat_names <- colnames(merged[5:22])
colnames(ibaq_nums)[5:22] <-goat_names


# get absolute numbers by dividing ibaq numbers by total protein loaded on column (measured by AAA)


ibaq.pd.nums.sums<- colSums(ibaq_nums[,5:22], na.rm=TRUE)
ibaq.pd.nums.factors <- 3/ibaq.pd.nums.sums
pd.ibaq.abs1 <-ibaq_nums[,5:22]
pd.ibaq.abs <-data.frame(mapply(`*`,pd.ibaq.abs1,ibaq.pd.nums.factors))
pd.ibaq.abs.final <-cbind(ibaq_nums[,1:4],pd.ibaq.abs)

# rearrage columns


# pd.ibaq.abs.final <-pd.ibaq.abs.final[,c(1:4,7,5,6,10,9,8,12,13,11,14,16,15,18,17,19,22,21,20)]

write.csv(pd.ibaq.abs.final,file="xtandem_ibaq_p_phos.tic.csv")

# load IBAQ data from maxquant and select goat proteins Need maxquant data for this

#mac WD
setwd("/Volumes/Lab$/service/on_campus/DePeters/Final Analysis/ibaq/combined/txt")

# PC WD
setwd("P:/service/on_campus/DePeters/Final Analysis/ibaq/combined/txt")

maxquant.raw <- read.csv("proteinGroups.csv",header=TRUE)
maxquant.raw.colnames <- c("Goat.5.1","Goat.5.3","Goat.6.3","Goat.1.2","Goat.2.2","Goat.3.1","Goat.3.2","Goat.1.3","Goat.2.1","Goat.4.1","Goat.6.2","Goat.4.3","Goat.1.1","Goat.4.2","Goat.2.3","Goat.3.3","Goat.5.2","Goat.6.1")

maxquant.ibaqs <-maxquant.raw[,c(1:10,146:163)]
colnames(maxquant.ibaqs)[11:28] <-maxquant.raw.colnames
write.csv(maxquant.ibaqs,file="maxquant.ibaqs")
maxquant.lfqs <-maxquant.raw[,c(1:10,164:181)]
colnames(maxquant.lfqs)[11:28] <-maxquant.raw.colnames
write.csv(maxquant.lfqs,file="maxquant.lfqs")

# sum maxquant ibaqs per column
maxquant.goat.ibaqs.sum<- colSums(maxquant.ibaqs[,11:28], na.rm=TRUE)
maxquant.goat.ibaqs.factors <- 3/maxquant.goat.ibaqs.sum

# divide ibaq numbers by factros 
maxquant.ibaq.abs1 <-maxquant.ibaqs[,11:28]
maxquant.ibaq.abs <-data.frame(mapply(`*`,maxquant.ibaq.abs1,maxquant.goat.ibaqs.factors))
maxquant.ibaq.abs.final <-cbind(maxquant.ibaqs[,1:10],maxquant.ibaq.abs)
maxquant.ibaq.abs.ordered.final <-maxquant.ibaq.abs.final[,c(1:10,23,14,18,19,15,25,16,17,26,20,24,22,11,27,12,28,21,13)]
write.csv(maxquant.ibaq.abs.ordered.final,file="maxquant_15_ibaq_abs_v2.csv")
