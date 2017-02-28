# plot bar plot of alignment statistic summary
# Rscript ~/cmd/alignplot.r bowtie2

args <- commandArgs(TRUE)
library(RColorBrewer)

alignplot=function(type)
{
if (type=="bowtie2")
{
fls=Sys.glob('*.log')
fls
dals=lapply(fls,read.delim,header=F,strip.white=T,stringsAsFactors=F)
dals
da=do.call(cbind,dals)
colnames(da)=nakename(fls)
da
da=da[6:9,]
da
rownames(da)=c('total','failed','unique','multiple')

tda=apply(da,1,sub,pattern=').*',replacement=')')
tda=cbind(sample=rownames(tda),tda)
tda
write.cn(tda,'align.summary.tsv')

tda=apply(tda[,-1],2,sub,pattern=' \\(.*',replacement='')
class(tda)='numeric'
tda
tda=cbind(tda,aligned=tda[,1]-tda[,2])

}
if (type=="tophat2")
{
fls=Sys.glob('*/align_summary.txt')
fls
dals=lapply(fls,read.delim,header=F,strip.white=T,stringsAsFactors=F)
dals
da=do.call(cbind,dals)
colnames(da)=nakename(dirname(fls))
da
da=da[2:4,]
da
rownames(da)=c('total','aligned','multiple')

tda=apply(da,1,gsub,pattern='.*: +| have.*| of input',replacement='')
tda=cbind(sample=rownames(tda),tda)
tda
write.cn(tda,'align.summary.tsv')

tda=apply(tda[,-1],2,sub,pattern=' \\(.*',replacement='')
class(tda)='numeric'
tda
tda=cbind(tda,failed=tda[,1]-tda[,2],unique=tda[,2]-tda[,3])

}

tmp=melt(tda)
head(tmp)
colnames(tmp)=c('sample','alignment','reads')
tmp$alignment=factor(tmp$alignment,levels=c('total','aligned','unique','multiple','failed'))
g=ggplot(tmp,aes(sample,reads,fill=alignment))+geom_bar(stat='identity',position='dodge')+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
g
ggsave('align.summary.png',g)
}

alignplot(type=args)