# takes scaffold sample report showing Total TIC
#calcualtes # IBAQ peptides from exported subset.fasta
# calcuales a IBAQ based on TIC
# calcualtes a percentatge of total sample (100% total)
# optional Calcualtes absolute amounts per peptide scalled to amount of protein loaded on column


#=================
# instructions
#=================
# Works with Total TIC or Total Precursor intensity
# Expot Subset Database Under Export drop-down
# for Total-TIC (if you don't have precursor intensity data).
#       Experiment -> Quant analysis -> quant method = Total TIC
#       Uncheck (use normalization)
#       Display Options -> Quantitative Value
#       Export -> Samples report -> Samples Report with clusters


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
setwd("/Users/bphinney/Desktop/")
#setwd("P:/service/on_campus/DePeters/")

# read scaffold samples report Change the filename here
infile <- "Samples Report With Clusters for fayBestDB.xls"

temp<-readLines(infile)
temp <- temp[which(temp!="")]
temp <- temp[which(temp!="END OF FILE")]
write.table(temp, file="tmpfile", row.names=FALSE, col.names=FALSE,
sep="\t", quote=FALSE)
tt <- readLines("tmpfile")
linetoskip <- grep("^#", tt) - 1
temp2 <- read.table(file="tmpfile", sep="\t", fill=TRUE, skip=linetoskip,
header=FALSE, comment.char="@", stringsAsFactors=FALSE)
rownames(temp2) <- temp2$V1
dat <- temp2[2:dim(temp2)[1],]
colnames(dat) <- temp2[1,]
colnames(dat) <-make.names(colnames(dat))
names(dat)[names(dat)=="X."] <- "Scaffold.Protein.Number"


infasta <- "Subset Database.fasta"


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
write.csv(ibaqs_all, file="Number of ibaq peptides per protein.csv")

dat$Accession.Number <- sub(" .*","",dat$Accession.Number)
#merge df's
merged<- merge(ibaqs_all,dat,by="Accession.Number")
merged$Protein.Grouping.Ambiguity <-NULL
merged$Quantitative.Variance <-NULL
merged$Taxonomy <-NULL
merged$Visible. <-NULL
merged$Starred. <-NULL
merged$NA. <-NULL
ibaq_nums<-merged[,1:5]

for (i in 6:ncol(merged)) {
    merged[,i] <-as.numeric(merged[,i])
}
#divide total peak are by IBAQ numbers
for (i in 1:nrow(merged)) {
    for (j in 6:ncol(merged))
    ibaq_nums[i,j] <- merged[i,j]/merged[i,2]
}

col_names <- colnames(merged[6:ncol(merged)])
colnames(ibaq_nums)[6:ncol(merged)] <-col_names
write.csv(ibaq_nums,file="Protein Total TIC divided by ibaq numbers.csv")

#remove reverses and trypsin and cluster members
revs<- grep("^REVERSE",ibaq_nums$Accession.Number)
ibaq_nums <- ibaq_nums[-revs,]
trypsin <- grep("TRYP",ibaq_nums$Accession.Number)
ibaq_nums <- ibaq_nums[-trypsin,]
#lysC <- grep("LYSC",ibaq_nums$Accession.Number)
#ibaq_nums <- ibaq_nums[-lysC,]
cluster.members<- grep("\\.",ibaq_nums$Scaffold.Protein.Number)
ibaq_nums <- ibaq_nums[-cluster.members,]
# calculate Percent of total protein in sample this is good if you don't knwo the total amoount loaded on the column

row_names<-ibaq_nums[,4]
ibaq_nums.ma <-as.matrix(ibaq_nums[,c(6:ncol(merged))])
colnames(ibaq_nums.ma) <-col_names
rownames(ibaq_nums.ma) <-row_names
percent_total <- prop.table(ibaq_nums.ma,2)
write.csv(percent_total,file="Protein Total TIC IBAQ Percent total sample clusters faybestDB.csv")


# get absolute numbers by dividing ibaq numbers by total protein loaded on column (measured by AAA). Only do this if you know how much you loaded on the column


#ibaq.pd.nums.sums<- colSums(ibaq_nums[,5:22], na.rm=TRUE)
#ibaq.pd.nums.factors <- 3/ibaq.pd.nums.sums
#pd.ibaq.abs1 <-ibaq_nums[,5:22]
#pd.ibaq.abs <-data.frame(mapply(`*`,pd.ibaq.abs1,ibaq.pd.nums.factors))
#pd.ibaq.abs.final <-cbind(ibaq_nums[,1:4],pd.ibaq.abs)

#write.csv(pd.ibaq.abs.final,file="pd14_ibaq_p.csv")

# load IBAQ data from maxquant and select goat proteins Need maxquant data for this

#mac WD
#setwd("/Volumes/Lab$/service/on_campus/DePeters/Final Analysis/ibaq/combined/txt")

# PC WD
#setwd("P:/service/on_campus/DePeters/Final Analysis/ibaq/combined/txt")

#maxquant.raw <- read.csv("proteinGroups.csv",header=TRUE)
#maxquant.raw.colnames <- c("Goat.5.1","Goat.5.3","Goat.6.3","Goat.1.2","Goat.2.2","Goat.3.1","Goat.3.2","Goat.1.3","Goat.2.1","Goat.4.1","Goat.6.2","Goat.4.3","Goat.1.1","Goat.4.2","Goat.2.3","Goat.3.3","Goat.5.2","Goat.6.1")

#maxquant.ibaqs <-maxquant.raw[,c(1:10,146:163)]
#colnames(maxquant.ibaqs)[11:28] <-maxquant.raw.colnames
#write.csv(maxquant.ibaqs,file="maxquant.ibaqs")
#maxquant.lfqs <-maxquant.raw[,c(1:10,164:181)]
#colnames(maxquant.lfqs)[11:28] <-maxquant.raw.colnames
#write.csv(maxquant.lfqs,file="maxquant.lfqs")

# sum maxquant ibaqs per column
#maxquant.goat.ibaqs.sum<- colSums(maxquant.ibaqs[,11:28], na.rm=TRUE)
#maxquant.goat.ibaqs.factors <- 3/maxquant.goat.ibaqs.sum

# divide ibaq numbers by factros
#maxquant.ibaq.abs1 <-maxquant.ibaqs[,11:28]
#maxquant.ibaq.abs <-data.frame(mapply(`*`,maxquant.ibaq.abs1,maxquant.goat.ibaqs.factors))
#maxquant.ibaq.abs.final <-cbind(maxquant.ibaqs[,1:10],maxquant.ibaq.abs)
#maxquant.ibaq.abs.ordered.final <-maxquant.ibaq.abs.final[,c(1:10,23,14,18,19,15,25,16,17,26,20,24,22,11,27,12,28,21,13)]
#write.csv(maxquant.ibaq.abs.ordered.final,file="maxquant_15_ibaq_abs_v2.csv")