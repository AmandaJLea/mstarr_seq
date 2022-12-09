############
# Merging counts from 24 samples
############

library(data.table)
lids<-read.delim('LID_info.txt')
b<-fread('/Genomics/ayroleslab2/alea/ref_genomes/hg38/hg38_200bp_windows.bed',header=F)

all_counts<-matrix(ncol=24,nrow=dim(b)[1])

for (p in 1:length(lids$LID)) {
  myfile<-paste0(lids$LID[p],'_200kb_cov.bed')
  a<-fread(myfile,header=F)
 identical(a$V2,b$V2);identical(a$V1,b$V1)
 all_counts[,p]<-a$V4
}

all_counts<-as.data.frame(all_counts)
names(all_counts)<-lids$LID
all_counts$site<-paste(b$V1,b$V2,sep='_')

############
# Filtering counts
############

# missing data
all_counts$zero_counts_DNA_un<-apply(all_counts[,subset(lids,type=='DNA' & condition=='unmeth')$LID],1,function(x) length(which(x==0)))
all_counts$zero_counts_RNA_un<-apply(all_counts[,subset(lids,type=='RNA' & condition=='unmeth')$LID],1,function(x) length(which(x==0)))

all_counts$zero_counts_DNA_m<-apply(all_counts[,subset(lids,type=='DNA' & condition=='meth')$LID],1,function(x) length(which(x==0)))
all_counts$zero_counts_RNA_m<-apply(all_counts[,subset(lids,type=='RNA' & condition=='meth')$LID],1,function(x) length(which(x==0)))

# filter for sites with DNA for at least 3 samples
all_counts2<-subset(all_counts,zero_counts_DNA_un<4 & zero_counts_DNA_m<4 & (zero_counts_RNA_m<4 | zero_counts_RNA_un<4))

# median coverage DNA
all_counts2$median_DNA_un<-apply(all_counts2[,subset(lids,type=='DNA' & condition=='unmeth')$LID],1,function(x) median(x[which(x>0)]))
all_counts2$median_DNA_m<-apply(all_counts2[,subset(lids,type=='DNA' & condition=='meth')$LID],1,function(x) median(x[which(x>0)]))

all_counts3<-subset(all_counts2, median_DNA_un>2 & median_DNA_m>2 & median_DNA_un<1001 & median_DNA_m<1001  )

write.table(all_counts3[,1:25],'7Dec22_mSTARR_counts_RAJ_K562.txt',row.names=F,sep='\t',quote=F)