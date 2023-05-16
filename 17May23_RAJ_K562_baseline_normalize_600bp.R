#Prepare voom normalized counts file
require(limma); require(edgeR); require(statmod); library(qvalue)
tt<-'nullwol31250' #IFN or dex or nullwol31250
zz<-'null' #IFN or dex or null
all<-read.table(paste0('rnadnacounts_3rna_6dna_highdnarepeat_dups_600bpwin_',tt,'.txt',sep=''),sep='\t')
info<-read.delim('~/Desktop/github_mstarr/1_mSTARR-seq_captures_enhancer_activity/coverage_report_Novaseq.txt')
info<-subset(info,libraryid!='L31250' & libraryid!='L31286')
info<-subset(info,treatment==zz) #Keep treatment of interest
info<-droplevels(info)
info$sampletype<-as.factor(info$sampletype)
all<-all[,info$libraryid]
info1<-info$meth
info2<-info$sampletype
design1<-model.matrix(~info2)
design1
all[is.na(all)]<-0
dge <- DGEList(counts=(all))
dge1 <- calcNormFactors(dge)
v <- voomWithQualityWeights(dge1,design=design1,plot=T)
nn<-v$E
write.table(nn,file='rnadnavoom_3rna_6dna_highdnarepeat_dups_600bpwin_nullwol31250.txt',sep='\t',quote=F)
