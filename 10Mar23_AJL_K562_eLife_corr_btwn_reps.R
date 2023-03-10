library(data.table)
data=fread('7Dec22_mSTARR_normcounts_AJL_K562_eLife.txt')

library(reshape2)
info<-as.data.frame(names(data))
names(info)<-'sample'
info$treatment<-c(rep('UN',12),rep('M',10))
info$type<-c(rep('DNA',6),rep('RNA',6),rep('DNA',5),rep('RNA',5))

# correlations between samples
melted_cormat <- melt(cor(data))
melted_cormat2 <- merge(melted_cormat,info,by.x='Var1',by.y='sample')
melted_cormat3 <- merge(melted_cormat2,info,by.x='Var2',by.y='sample')

# high correlations between DNA replicates
library(ggplot2)
dna<-subset(melted_cormat3,type.x=='DNA' & type.y=='DNA' & value<1)
ggplot(data = dna, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+theme_bw(13)+scale_fill_gradient(low = "yellow", high = "blue", na.value = NA)

# lower correlations between RNA replicates (expected); especially methylated replicates
rna<-subset(melted_cormat3,type.x=='RNA' & type.y=='RNA' & value<1)
ggplot(data = rna, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+theme_bw(13)+scale_fill_gradient(low = "yellow", high = "blue", na.value = NA)

# low correlations among RNA reps are driven by low correlations when there are low counts/weak enhancers
smoothScatter(data$UN3.y,data$UN4.y)
  
# mean counts aggregated across replicates are still consistently different between meth and unmeth for MD enhancers
# the noisiness does mean that when we call an MD enhancer, it's something that's "on" in the unmethylated condition and usually almost or completely "off" in the methylated condition
res=fread('7Dec22_results_AJL_K562_eLife.txt')
enh<-subset(res,info.y=='MD')
smoothScatter(enh$norm_RNA_unmeth,enh$norm_RNA_meth)
x=c(0,1);y=c(0,1);abline(lm(y~x),lty=2)

