library("GenomicRanges")
library("GenomicAlignments")
library(parallel)

####################################################################################
###################################### Peaks #######################################
peak.gr=function(peak)
{
if (colnames(peak)[1]== "chr")
peakgr=GRanges(seqnames=peak$chr,ranges=IRanges(start=peak$start,end=peak$end),strand=rep("*",nrow(peak)),peak_index=rownames(peak)) else
if (colnames(peak)[1]== "seqnames")
peakgr=GRanges(seqnames=peak$seqnames,ranges=IRanges(start=peak$peak_start,end=peak$peak_end),strand=rep("*",nrow(peak)),peak_index=peak$peak_index) else
peakgr=GRanges(seqnames=peak[,1],ranges=IRanges(start=peak[,2],end=peak[,3]),strand=rep("*",nrow(peak)),peak_index=rownames(peak))
if ("summit" %in% names(peak)) peakgr$summit=peak$summit
peakgr
}

olmerge=function(peak.gr1,peak.gr2)
{
    ol1=subsetByOverlaps( peak.gr1,peak.gr2, maxgap = 0L, minoverlap = 1L, type ="any", ignore.strand = T)
    ol2=subsetByOverlaps( peak.gr2,peak.gr1, maxgap = 0L, minoverlap = 1L, type ="any", ignore.strand = T)
    reduce(c(ol1,ol2))
}

olvenn=function(pgrls,maxgap=0,figure=F)
{
	
	x=GRangesList(pgrls)
	p.merge=reduce(unlist(x))
	
	f=function(pgr,p.merge)
	{
	ol=findOverlaps( p.merge, pgr, maxgap = maxgap, minoverlap = 1L, type ="any", select ="all", ignore.strand = T)
	(1:length(p.merge)) %in% queryHits(ol)
	}
	
	b=lapply(pgrls,f,p.merge)
	b.mx=do.call(cbind,b)
	#colnames(b.mx)=names(pgrls)
	n=1:length(p.merge)
	mcols(p.merge)=cbind(as.data.frame(mcols(p.merge)),as.data.frame(n),b.mx)
	
	if (figure==T)
	{
	library(limma)
    ##vennCounts(b.mx)
	x=vennCounts(as.data.frame(mcols(p.merge)[,-1]))
	vennDiagram(x)
	if (length(pgrls)==2)
		{
		myvenn(x[4,3]+x[3,3],x[4,3]+x[2,3],x[4,3],names(pgrls),paste(c("venn",names(pgrls),"png"),collapse='.')) ## if only two peaks overlap
		}
	}
	
	p.merge
}

# x=vennCounts(as.data.frame(mcols(pmerge.tss)[,2:3]))
#venn(x[4,3]+x[3,3],x[4,3]+x[2,3],x[4,3],c("H3K4me3","MLL2"),"venn.k4tss.mll2tss.png")
myvenn=function(a,b,ab,name=c("A","B"),output,imagetype="png") # only "png","tiff","svg"
{
  require(VennDiagram)
### venn.diagram(list(B = 1:1800, A = 1571:2020),fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =2, fontfamily =3, filename = "trial2.emf")
   cls=list(1:a, (a-ab+1):(a+b-ab))
   names(cls)=name
   venn.diagram(cls,fill = c("red", "green"), alpha = c(0.5, 0.5), cex = 2,cat.fontface = 4,lty =1, fontfamily =3, filename = output, imagetype=imagetype, margin=0.05)
} 

#hypergeometric test for venn gene lists
venn.test=function(a,b,ab,n)
{
phyper(ab,a,n-a,b,lower.tail=F)

  #length(unique(mm9ref$Ensembl.Gene.ID)) #20130
   # two gene list
   #m=805
   #k=1502
   # intersect
   #q=497
   #n=20130-m
   #phyper(q,m,n,k, lower.tail=F) # 0
}

#### annpeak
annpeak=function(peak.gr,tr.gr,tssup=0,tssdown=0)
{
library(GenomicRanges)

#olann=distanceToNearest(peak.gr,tr.gr,select="all",ignore.strand = TRUE)
olann=distanceToNearest(peak.gr,tr.gr,ignore.strand = TRUE)
olann=as.data.frame(olann)
ann=cbind(as.data.frame(peak.gr),as.data.frame(tr.gr)[olann$subjectHits,],distance=olann$distance) #bioconductor <2.12
names(ann)[2]="peak_start"
names(ann)[3]="peak_end"
names(ann)[4]="peak_width"
names(ann)[5]="peak_strand"

ann$peak_to_feature="body"
up=which((ann$strand=="+" & ann$start > ann$peak_end  ) |(ann$strand=="-" & ann$peak_start > ann$end ))
if (length(up) != 0) ann[up,]$peak_to_feature="upstream"
down=which((ann$strand=="+" & ann$peak_start > ann$end ) |(ann$strand=="-" & ann$start > ann$peak_end ))
if (length(down) != 0) ann[down,]$peak_to_feature="downstream"

## correct tss overlap
s<- as(strand(tr.gr),"vector")=="+"
tss=ifelse(s,start(tr.gr),end(tr.gr))
tss.gr=GRanges(seqnames=seqnames(tr.gr),ranges=IRanges(start=ifelse(s,tss-tssup,tss-tssdown),end=ifelse(s,tss+tssdown,tss+tssup)),strand(tr.gr),tr.id=tr.gr$tr.id)

olann_tss=distanceToNearest(peak.gr,tss.gr,ignore.strand = TRUE)
olann_tss=as.data.frame(olann_tss)
ann_tss=cbind(as.data.frame(peak.gr),as.data.frame(tr.gr)[olann_tss$subjectHits,],distance=olann_tss$distance)
ann_tss$peak_to_feature="overlap_tss"

oltss=which(ann_tss$distance==0)
if (length(oltss)!=0) ann[oltss,]=ann_tss[oltss,]

ann$tss=ifelse(ann$strand=="+",ann$start,ann$end)
#ann$distance_summit_to_tss=ifelse(ann$strand=="+", (ann$summit+ann$peak_start)-ann$tss, ann$tss - (ann$summit+ann$peak_start))	#use for sorting
#ann$distance_to_gene=ifelse(ann$peak_to_feature=="upstream", -ann$distance, ann$distance)
ann$center=round((ann$peak_start+ann$peak_end)/2)
ann$distance_center_to_tss=ifelse(ann$strand=="+", ann$center-ann$tss, ann$tss - ann$center)	#use for sorting
#tmp=which(ann$peak_to_feature=="overlap_tss")
#if (length(tmp)!=0) 
#{ann[tmp,]$distance_summit_to_tss= 0
#ann[tmp,]$distance_center_to_tss= 0
#}

ann

}

plotpie=function(annp,main,out.pdf=F)
{
mycolor <- c("#7fc1f9","#fe9e6e","#b5ea7c","#d682fd","#676667")
if (out.pdf==T) pdf(paste0("pie.",main,".pdf")) else
png(paste0("pie.",main,".png"))
par(mar=c(0,0,4,0))
n= summary(as.factor(annp$peak_to_feature))[c(3,1,2,4)] 
pie(n,labels="",col=mycolor,border=NA,clockwise = T,init.angle=180,main=main)
points(0,0,pch=16,cex=par("cex")*35,col="white")
N=nrow(annp)
l=paste0(names(n)," = ",round(n/N*100),"%")
l=c(paste0("N = ",N),"peak to features:",l)
legend(-0.4,0.3,legend=l,ncol=1,col=c("white","white",mycolor),pch=15,bty="n",cex=1.3)
dev.off()
}

peakpie <- function(annp,fname)
{
library(ggplot2)
annp$peak_to_feature=factor(annp$peak_to_feature,levels=levels(as.factor(annp$peak_to_feature))[c(3,1,2,4)])
ggplot(annp,aes(1,fill=peak_to_feature))+geom_bar()+coord_polar(theta="y")
}

annpeak_bytss=function(peak.gr,tr.gr,tssup=0,tssdown=0)
{
library(GenomicRanges)

s<- as(strand(tr.gr),"vector")=="+"
tss=ifelse(s,start(tr.gr),end(tr.gr))
tss.gr=GRanges(seqnames=seqnames(tr.gr),ranges=IRanges(start=ifelse(s,tss-tssup,tss-tssdown),end=ifelse(s,tss+tssdown,tss+tssup)),strand(tr.gr),tr.id=tr.gr$tr.id)

#olann_tss=distanceToNearest(peak.gr,tss.gr,select="all",ignore.strand = TRUE)
olann_tss=distanceToNearest(peak.gr,tss.gr,ignore.strand = TRUE)
olann_tss=as.data.frame(olann_tss)
ann_tss=cbind(as.data.frame(peak.gr),as.data.frame(tr.gr)[olann_tss$subjectHits,],distance=olann_tss$distance)
names(ann_tss)[2]="peak_start"
names(ann_tss)[3]="peak_end"
names(ann_tss)[4]="peak_width"
names(ann_tss)[5]="peak_strand"

ann_tss$tss=ifelse(ann_tss$strand=="+",ann_tss$start,ann_tss$end)
#ann_tss$distance_summit_to_tss=ifelse(ann_tss$strand=="+", (ann_tss$summit+ann_tss$peak_start)-ann_tss$tss, ann_tss$tss - (ann_tss$summit+ann_tss$peak_start))	#use for sorting
#ann_tss[which(ann_tss$distance==0),]$distance_summit_to_tss= 0
ann_tss$center=round((ann_tss$peak_start+ann_tss$peak_end)/2)
ann_tss$distance_center_to_tss=ifelse(ann_tss$strand=="+", ann_tss$center-ann_tss$tss, ann_tss$tss - ann_tss$center)
if (min(ann_tss$distance)== 0) ann_tss[which(ann_tss$distance==0),]$distance_center_to_tss= 0	
ann_tss

}

plotpie_tss=function(annp_bytss,main,out.pdf=F)
{
if (out.pdf==T) pdf(paste0("pie.",main,".pdf")) else
png(paste0("pie.",main,".png"))
par(mar=c(0,0,4,0))
  x=annp_bytss$distance_center_to_tss/1000
  g=cut(x,breaks=c(min(x),-50,-10,-5,-1,0,1,5,10,50,max(x)))
  xg=split(x,g)
  n=unlist(lapply(xg,length))
  dis=n/sum(n)  
  dis2=dis[6:10]+dis[5:1]
  pie(dis2,labels="",col=rev(brewer.pal(5,"Blues")),border=NA,clockwise = T,init.angle=180,main=main)
  points(0,0,pch=16,cex=par("cex")*35,col="white")
  #legend("topright",legend=c("<1 kb","1-5 kb","5-10 kb","10-50 kb",">50 kb"),pch=15,cex=1.5,col=rev(brewer.pal(5,"Blues")))
  N=nrow(annp_bytss)
l=paste0(c("<1 kb","1-5 kb","5-10 kb","10-50 kb",">50 kb")," : ",round(dis2*100),"%")
l=c(paste0("N = ",N),"peak to tss:",l)
legend(-0.3,0.4,legend=l,ncol=1,col=c("white","white",rev(brewer.pal(5,"Blues"))),pch=15,bty="n",cex=1.3)
dev.off()
  
 }
 
plotbar=function(annp_bytss.ls)
{
  f=function(annp_bytss)
 {
  x=annp_bytss$distance_center_to_tss/1000
  g=cut(x,breaks=c(min(x),-50,-10,-5,-1,0,1,5,10,50,max(x)))
  xg=split(x,g)
  n=unlist(lapply(xg,length))
  n/sum(n)
  }
  
  if(length(annp_bytss.ls)==1) dis=f(annp_bytss.ls[[1]]) else
  dis=do.call(cbind,lapply(annp_bytss.ls,f))

  barplot(as.matrix(dis),col=c(brewer.pal(5,"Blues"),rev(brewer.pal(5,"Blues"))),horiz=T,xaxt="n",xlim=c(0,1.5))
  #axis(1,at=sum(dis[1:5]),labels="TSS")
  legend("topright",legend=c("<1 kb","1-5 kb","5-10 kb","10-50 kb",">50 kb"),pch=15,cex=1.5,col=rev(brewer.pal(5,"Blues")))
  
##     f=function(annp_bytss)
##  {
##   x=annp_bytss$distance_center_to_tss/1000
##   g=cut(x,breaks=c(min(x),-50,-10,-5,-1,0,1,5,10,50,max(x)))
##   xg=split(x,g)
##   names(xg)[1]="<= -50"
##   names(xg)[10]=">= 50"
##   tmp=melt(xg)
##   tmp$L1
##   }
##   
##   dis=lapply(annp.tss,f)
##   dis.m=melt(dis)
##   #dis.m$name=dis.m$value
##   #levels(dis.m$name)=c("<1 kb","<1 kb","5-10 kb","10-50 kb","1-5 kb",">50 kb",">50 kb","10-50 kb","1-5 kb","5-10 kb")
##   dis.m$value=factor(dis.m$value,levels=levels(dis.m$value)[c(6,8,3,9,2,1,5,10,4,7)])
##   
##   ggplot(dis.m,aes(L1,fill=value)) + 
##   geom_bar(position="fill",color="black") +
##   scale_fill_manual(values=c(brewer.pal(5,"Blues"),rev(brewer.pal(5,"Blues"))))+
##   coord_flip()+
##   xlab("") +
##   ylab("") +  
##   scale_y_continuous(breaks=c(0,0.5,1), labels=c("upsream","TSS","downstream"))
##   
 }
 
# anngene by peak
anngene=function(peak.gr,tr.gr)
{
ann1k=annpeak(peak.gr,tr.gr,tssup=1000,tssdown=1000)
ann5k=annpeak(peak.gr,tr.gr,tssup=5000,tssdown=5000)
ann10k=annpeak(peak.gr,tr.gr,tssup=10000,tssdown=10000)

gene.tss1k=unique(ann1k$gene.id[ann1k$peak_to_feature=="overlap_tss"])
gene.tss5k=unique(ann5k$gene.id[ann5k$peak_to_feature=="overlap_tss"])
gene.tss10k=unique(ann10k$gene.id[ann10k$peak_to_feature=="overlap_tss"])
gene.body=unique(ann1k$gene.id[ann1k$peak_to_feature=="body"])

tr.gr$peakatTss=">10kb"
tmp=which(tr.gr$gene.id %in% gene.tss10k) ## asign order is important
tr.gr$peakatTss="<=10kb"
tmp=which(tr.gr$gene.id %in% gene.tss5k)
tr.gr$peakatTss="<=5kb"
tmp=which(tr.gr$gene.id %in% gene.tss1k)
tr.gr$peakatTss="<=1kb"
tr.gr$peakatBody=tr.gr$gene.id %in% gene.body

tr.gr
}

## gene nearest peak
nearest.peak=function(peak.gr,tr.gr,tssup=0,tssdown=0)
{
s<- as(strand(tr.gr),"vector")=="+"
tss=ifelse(s,start(tr.gr),end(tr.gr))
tss.gr=GRanges(seqnames=seqnames(tr.gr),ranges=IRanges(start=ifelse(s,tss-tssup,tss-tssdown),end=ifelse(s,tss+tssdown,tss+tssup)),strand(tr.gr),id=names(tr.gr))

olann=distanceToNearest(tss.gr,peak.gr,ignore.strand = TRUE)
olann=as.data.frame(olann)
ann=cbind(as.data.frame(tss.gr)[olann$queryHits,],as.data.frame(peak.gr)[olann$subjectHits,],distance=olann$distance) #bioconductor <2.12
colnames(ann)[7:11]=paste0('peak_',colnames(ann)[7:11])
ann
}




