############
# All files are from the eLife 2018 dataset of K562s assayed at baseline/unstimulated
############

# module load R/4.2.0

# all counts from less conservative pipeline
library(data.table)
setwd("/Genomics/ayroleslab2/alea/mSTARR/mSTARR_AI")
all2=fread('7Dec22_mSTARR_counts_AJL_K562.txt')

# eLife regions
enh1=read.delim('Supplementary_File_3_eLife.txt')
MD1=read.delim('Supplementary_File_4_eLife.txt')
enh1$site<-paste(enh1[,7],enh1[,8],sep='_')
MD1$site<-paste(MD1[,10],MD1[,11],sep='_')

all3=subset(all2,site %in% enh1$site)

############
# Normalize counts
############

library(DESeq2); library(edgeR); library(qvalue)
treatment<-c(rep('1_unmeth',12),rep('2_meth',10))
type<-c(rep('1_DNA',6),rep('2_RNA',6),rep('1_DNA',5),rep('2_RNA',5))
design<-model.matrix(~treatment + treatment:type)
dge <- DGEList(counts=(all3[,c(2:23),with=F]))
dge1 <- calcNormFactors(dge)
v <- voomWithQualityWeights(dge1,design=design,plot=FALSE)
write.table(v$E,'7Dec22_mSTARR_normcounts_AJL_K562_eLife.txt',row.names=F,sep='\t',quote=F)

############
# Summary of normalized counts
############

norm=fread('7Dec22_mSTARR_normcounts_AJL_K562_eLife.txt')
norm_DNA_unmeth<-apply(norm[,1:6],1,mean)
norm_RNA_unmeth<-apply(norm[,7:12],1,mean)
norm_DNA_meth<-apply(norm[,13:17],1,mean)
norm_RNA_meth<-apply(norm[,18:22],1,mean)
norm_all<-as.data.frame(cbind(norm_DNA_unmeth,norm_RNA_unmeth,norm_DNA_meth,norm_RNA_meth))
norm_all$id<-1:dim(norm_all)[1]
write.table(norm_all,'7Dec22_mSTARR_normmeans_AJL_K562_eLife.txt',row.names=F,sep='\t',quote=F)

############
# Identify regulatory activity - RNA>DNA in at least 1 condition 
############

design<-model.matrix(~treatment + treatment:type)
fit <-lmFit(v,design)
fit_dna <- eBayes(fit)
fit_dna <- as.data.frame(cbind(fit_dna$p.value[,3:4],fit_dna$coefficient[,3:4]))
fit_dna$id<-1:dim(norm_all)[1]
names(fit_dna)<-c('pval_unmeth','pval_meth','beta_unmeth','beta_meth','id')
fit_dna$q_unmeth<-qvalue(fit_dna[,1])$qvalues
fit_dna$q_meth<-qvalue(fit_dna[,2])$qvalues
fit_dna_keep<-subset(fit_dna,(q_unmeth<0.1 & beta_unmeth>0 ) |(q_meth<0.1 & beta_meth>0 )); dim(fit_dna_keep)
fit_dna$info<-'not_regulatory'
fit_dna$info[which(fit_dna$id %in% fit_dna_keep$id)]<-'regulatory'
write.table(fit_dna,'7Dec22_model1_AJL_K562_eLife.txt',row.names=F,sep='\t',quote=F)

############
# Identify MD regulatory activity - RNA/DNA x unmeth/meth effect 
############

design<-model.matrix(~treatment*type)
fit <-lmFit(v,design)
fit_int <- eBayes(fit)
fit_int <- as.data.frame(cbind(fit_int$p.value[,2:4],fit_int$coefficient[,2:4]))
names(fit_int)<-c('pval_treat','pval_type','pval_int','beta_treat','beta_type','beta_int')
fit_int$id<-1:dim(fit_int)[1]
fit_int2<-subset(fit_int,id %in% fit_dna_keep$id)
fit_int2$int_qval<-qvalue(fit_int2[,3])$qvalues
length(which(fit_int2$int_qval<0.1))
length(which(fit_int2$int_qval<0.1 & fit_int2[,6]<0))
fit_int2$info<-'not_MD'
fit_int2$info[which(fit_int2$int_qval<0.1)]<-'MD'
write.table(fit_int2,'7Dec22_model2_AJL_K562_eLife.txt',row.names=F,sep='\t',quote=F)

############
# Combine and check files
############

library(data.table)

# all counts from less conservative pipeline
library(data.table)
setwd("/Genomics/ayroleslab2/alea/mSTARR/mSTARR_AI")
all2=fread('7Dec22_mSTARR_counts_AJL_K562.txt')

# eLife regions
enh1=read.delim('Supplementary_File_3_eLife.txt')
MD1=read.delim('Supplementary_File_4_eLife.txt')
enh1$site<-paste(enh1[,7],enh1[,8],sep='_')
MD1$site<-paste(MD1[,10],MD1[,11],sep='_')

all3=subset(all2,site %in% enh1$site)
all3$id<-1:dim(all3)[1]
sites=as.data.frame(all3[,c('site','id'),with=F])

# new files
enh2=fread('7Dec22_model1_AJL_K562_eLife.txt')
MD2=fread('7Dec22_model2_AJL_K562_eLife.txt')
mean=fread('7Dec22_mSTARR_normmeans_AJL_K562_eLife.txt')

# merge new files
all=merge(enh2,MD2[,c('pval_int','beta_int','int_qval','info','id')],by='id',all.x=T)
all=merge(all,mean,by='id')
all$norm_ratio_unmeth<-all$norm_RNA_unmeth-all$norm_DNA_unmeth
all$norm_ratio_meth<-all$norm_RNA_meth-all$norm_DNA_meth
all=merge(all,sites,by='id')
all=as.data.frame(all[order(all$id),])
write.table(all,'7Dec22_results_AJL_K562_eLife.txt',row.names=F,sep='\t',quote=F)

# check against old files
enh1$site<-paste(enh1[,7],enh1[,8],sep='_')
MD1$site<-paste(MD1[,10],MD1[,11],sep='_')
tmp<-subset(enh1, (Effect.of.library.type.in.the.methylated.condition..FDR.corrected.p.value.<0.1 & Effect.of.library.type.in.the.methylated.condition..beta.<0 ) |(Effect.of.library.type.in.the.unmethylated.condition..FDR.corrected.p.value.<0.1 & Effect.of.library.type.in.the.unmethylated.condition..beta.<0 ) )

length(which(tmp$site %in% subset(all,info.x=='regulatory')$site))

tmp<-subset(MD1, (Effect.of.condition.x.library.type.interaction..FDR.corrected.p.value.<0.1 ) )

length(which(tmp$site %in% subset(all,info.y=='MD')$site))

