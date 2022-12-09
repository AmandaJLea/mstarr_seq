############
# All files are from the eLife 2018 dataset of K562s assayed at baseline/unstimulated
# Start with merging and filtering mSTARR-seq files
############

library(data.table)
meth_dna=fread('16Jan16_pileup_DNA_meth_v2.txt')
meth_rna=fread('16Jan16_pileup_RNA_meth_v2.txt')
unmeth_dna=fread('16Jan16_pileup_DNA_unmeth_v2.txt')
unmeth_rna=fread('16Jan16_pileup_RNA_unmeth_v2.txt')

unmeth_dna2<-unmeth_dna[,.(site,UN1,UN2,UN3,UN4,UN5,UN6)]
unmeth_rna2<-unmeth_rna[,.(site,UN1,UN2,UN3,UN4,UN5,UN6)]

meth_dna2<-meth_dna[,.(site,M1,M3,M4,M5,M6)]
meth_rna2<-meth_rna[,.(site,M1,M3,M4,M5,M6)]

all_unmeth<-as.data.frame(merge(unmeth_dna2,unmeth_rna2,by='site'))
all_meth<-as.data.frame(merge(meth_dna2,meth_rna2,by='site'))

# missing data
all_unmeth$zero_counts_DNA_un<-apply(unmeth_dna2[,.(UN1,UN2,UN3,UN4,UN5,UN6)],1,function(x) length(which(is.na(x))))
all_unmeth$zero_counts_RNA_un<-apply(unmeth_rna2[,.(UN1,UN2,UN3,UN4,UN5,UN6)],1,function(x) length(which(is.na(x))))
all_meth$zero_counts_DNA_m<-apply(meth_dna2[,.(M1,M3,M4,M5,M6)],1,function(x) length(which(is.na(x))))
all_meth$zero_counts_RNA_m<-apply(meth_rna2[,.(M1,M3,M4,M5,M6)],1,function(x) length(which(is.na(x))))

# filter for sites with DNA for at least 3 samples
all_unmeth2<-subset(all_unmeth,zero_counts_DNA_un<4 )
all_meth2<-subset(all_meth,zero_counts_DNA_m<3 )

# median coverage DNA
all_unmeth2$median_DNA_un<-apply(all_unmeth2[,c('UN1.x','UN2.x','UN3.x','UN4.x','UN5.x','UN6.x')],1,function(x) median(x,na.rm=TRUE))
all_meth2$median_DNA_m<-apply(all_meth2[,c('M1.x','M3.x','M4.x','M5.x','M6.x')],1,function(x) median(x,na.rm=TRUE))

all<-merge(all_unmeth2,all_meth2,by='site')
all2<-subset(all, median_DNA_un>2 & median_DNA_m>2 & median_DNA_un<1001 & median_DNA_m<1001 & (zero_counts_RNA_m<3 | zero_counts_RNA_un<4) )

treatment<-c(rep('unmeth',12),rep('meth',10))
type<-c(rep('DNA',6),rep('RNA',6),rep('DNA',5),rep('RNA',5))
rep<-c(1:6,1:6,7:11,7:11)

library('limma')
library('edgeR')
library('qvalue')

all2[is.na(all2)]<-0
write.table(all2[,c(1:13,19:28)],'7Dec22_mSTARR_counts_AJL_K562.txt',row.names=F,sep='\t',quote=F)