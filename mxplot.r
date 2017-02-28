
########################### calculate Heatmap matrix and plot metaplot and heatmap ##################
# query can be cover.rpm (for coverage) or peak.gr (for enrichment)
# mclapply(rpmls,mxplot,anntype="anncenter",ann.df=annmll_19_loss,left=-5000,right=5000,bin_width=25,mc.cores=6,mc.set.seed = FALSE)

mxplot=function(query,ann,anntype="tss",left=-5000,right=5000,bin_width=25)
{	
	library(GenomicRanges)
	library("GenomicAlignments")
	if (class(ann)=='data.frame') ann.gr=df2gr(ann,type=anntype) else ann.gr=ann
	
	#enrichment or coverage plot
	enrich <- class(query)=="GRanges"
	#if (class(query)=="GRanges") enrich=TRUE 
	if (enrich) cvg=coverage(query) else cvg=query	
	
	seql=unlist(lapply(cvg,length))
	seqlengths(ann.gr)=seql[seqlevels(ann.gr)]
	
	#resize window
	ann.resize.gr=promoters(ann.gr,upstream=abs(left),downstream=right)
	ann.gr=trim(ann.resize.gr)
	
	message("calculating")
	cvgls=cvg[ann.gr]
	strand=as.vector(strand(ann.gr))
	f=function(i)
	{
		vtmp=cvgls[[i]]
		n=length(right-left+1-length(vtmp))
		if (start(ann.gr)[i]<1) vtmp=append(vtmp,rep(0,n),after=0)
		
		if (end(ann.gr)[i]> seql[as.vector(seqnames(ann.gr[i]))] ) vtmp=append(vtmp,rep(0,n))
		
		if (strand[i] == -1 | strand[i] == '-') vtmp=rev(vtmp)
		
		if (enrich | bin_width==1) v=bin(as.vector(vtmp),bin_width) else v=bin_mean(as.vector(vtmp),bin_width)
		v
	}
	
	vls=lapply(1:length(ann.gr),f)
	v=do.call(rbind,vls)
	v
		
}
	
bin=function(v,bin)
{
	bincol=seq(1,length(v),bin)
	v[bincol]
}

bin_mean=function(v,bin)
{
	colMeans(array(v,c(bin,length(v)/bin)))
}

if (FALSE)
{#bin_mean for matirx
bin_mean=function(v,bin)
{
	bincol_start=seq(1,length(v),bin)
	bincol_end=bincol_start+bin-1
	v_bin=rep(0,length(v)/bin)
	for (i in 1:length(v_bin))
	{
		v_bin[i]=mean(v[bincol_start[i]:bincol_end[i]])
	}
	v_bin
}
}

reassign_log2=function(ratio)
{# reassign log2 ratio matrix for plot
tmp=ratio
tmp[tmp=="NaN"]=1
tmp[tmp== Inf]= max(tmp[tmp< Inf])
ratio=tmp
log2ratio=log2(ratio)
tmp=log2ratio
#tmp[tmp=="NaN"]=1
tmp[tmp== -Inf]= min(tmp[tmp> -Inf])
log2ratio=tmp
log2ratio
}

plotmx = function(mxls,output,cutoff=c(0,1),color=c("white","blue"),plot_win=c(0,20,70),win_label=c('-200','TSS','+500'),notes='')
{##### plot matrix in R
pal <- colorRampPalette(color)
x=do.call(cbind,mxls)
x=apply(x,2,rev)
x=t(x)

#increase size can make it smoother, save as pdf can make it like in java treeview
pdf(output,width=length(mxls)*2,height=10)
layout(matrix(c(rep(1,12),2,3), nrow = 7, ncol = 2, byrow = TRUE))

par(mar=c(1,1,10,1))
image(x, col=pal(250),breaks=c(min(x),seq(cutoff[1],cutoff[2],length=249),max(x)),useRaster=T,xaxt="n",yaxt="n")
#abline(v=1/length(mxls)*seq(0,length(mxls)))
abline(v=1/length(mxls)*seq(1,length(mxls)-1))

n=0.5/length(mxls)
text(seq(n,1-n,by=n*2),par("usr")[2]+0.1,labels=names(mxls),srt=30,pos=1,xpd=T,cex=2)

axis(side=1,at=1/nrow(x)*plot_win,labels=win_label,cex.axis=1.5 )


## color bar
mincol=color[1]
maxcol=color[length(color)]
col=pal(200)
par(mar=c(6,2,3,4))
if (cutoff[1]==0)
{
barplot(rep(1,300),col=c(col,rep(maxcol,100)),space=0,border=F,axes=F,xlab="reads per million (RPM)",cex.lab=2)
axis(side=1,at=c(0,50,100,150,200,300),labels=c(0,0.25,0.5,0.75,1,round(max(x))),cex.axis=1.5 )
}

if (cutoff[1]<0)
{
barplot(rep(1,400),col=c(rep(mincol,100),col,rep(maxcol,100)),space=0,border=F,axes=F,xlab="log2FC",cex.lab=1.5)
axis(side=1,at=c(0,100,200,300,400),labels=c(round(min(x)),cutoff[1],0,cutoff[2],round(max(x))),cex.axis=1.5 )
}

## text
par(mar=c(4,1,1,1))
plot.new()
text(0.5,0.2,labels=notes,cex=3)

dev.off()

}
 
metaplot=function(mxls,row.select,ylim=c(0,1),center="peak",color=NA,win="5 kb")
{
#select: select rows for calculate colMeans (used for grouping data)
ncols=ncol(mxls[[1]])
sub.mxls=lapply(mxls,"[",row.select,1:ncols)
meanls=lapply(sub.mxls,colMeans)
meanmx=do.call(cbind,meanls)

tmp=1:length(sub.mxls)
mycolor=ifelse(is.na(color),tmp,color)
#print(mycolor)

matplot(1:ncol(sub.mxls[[1]]),meanmx,pch=20,type="l",lty=1,lwd=2,col=mycolor,xlab="",ylab="", xaxt="n",ylim=ylim)
legend("topleft",names(sub.mxls),col=mycolor,lty=1)

if (center=="peak") axis(1, at=c(0,ncols/2,ncols),labels=c(paste0("-",win),"Peak center",win), col.axis="black", las=1) else
axis(1, at=c(0,ncols/2,ncols),labels=c(paste0("-",win),"TSS",win), col.axis="black", las=1)

#sd1=apply(sub.mxls[[1]],2,sd)
# x=1:400
# y1=meanls[[1]]+sd1
# y2=meanls[[1]]-sd1
# polygon(c(x,rev(x)),c(y2,rev(y1)),col="skyblue")
# lines(meanls[[1]],type="l")

}

if (FALSE)
{#### metaplot with confident interval
    f=function(v)
        {tmp=boxplot.stats(v)
         y=c(tmp$stats[2:4],tmp$conf)
         names(y)=c("25","median","75","confl","confu")
         y
         }

    x=metamx.t

    tmp=apply(rbind(x[[1]],x[[3]]),2,f)
    plot(tmp["median",],type="l",lwd=3,ylim=c(0,6))
    polygon(c(1:160,rev(1:160)),c(tmp["confu",],rev(tmp["confl",])),col="lightgrey",border=NA)
    #plot(tmp["confu",],type="l",lty=2)
    lines(tmp["median",],type="l",lwd=3)
    #lines(tmp["confl",],type="l",lty=2)

    tmp2=apply(rbind(x[[2]],x[[3]]),2,f)
    polygon(c(1:160,rev(1:160)),c(tmp2["confu",],rev(tmp2["confl",])),col="#FF000020",border=NA)
    #lines(tmp2["confu",],type="l",col="red",lty=2)
    lines(tmp2["median",],type="l",col="red",lwd=3)
    #lines(tmp2["confl",],type="l",col="red",lty=2)
}	

row_bin_mean=function(mx,bin) 
{#bin the rows to low down the size of matrix
frow=floor(nrow(mx)/bin) #throw away the last several lines
ind=gl(frow,bin)
#ind=as.factor(c(1,1,2,2,3))
mx.ls=split(mx[1:(frow*bin),],ind)
mx.reshape.ls=lapply(mx.ls,matrix,nrow=bin)
mx.ls2=lapply(mx.reshape.ls,colMeans)
do.call(rbind,mx.ls2)
}


########### old ############
##!!!!!! the following stop updating

########################### calculate Chip-seq occupancy from bigwig file  ##################

occupancy=function(anntype,ann.df,cover.rpm,left=-500,right=500)
{ ###!!!!!  use count.peak instead

library("GenomicRanges")

	#find type of ann.df
	switch(anntype,
		#use the full peak region
	annpeak={ann=data.frame(seqnames=ann.df$seqnames,start=ann.df$peak_start,end=ann.df$peak_end,strand=rep(1,nrow(ann.df)))},
	peak={ann=data.frame(seqnames=ann.df$chr,start=ann.df$start,end=ann.df$end,strand=rep(1,nrow(ann.df)))},
		#centred at tss, use mm9tss or hg19tss as ann.df
	ensembl={ann=data.frame(seqnames=ann.df$Chromosome.Name,start=ann.df$Transcript.Start..bp.,end=ann.df$Transcript.End..bp.,strand=ann.df$Strand)},
		#centred at tss, use annpeak output or as.data.frame(tr.gr) as ann.df
	tss={ann=ann.df},
		#centered at peak summit,use annpeak output as ann.df
	summit={ann=data.frame(seqnames=ann.df$seqnames,start=ann.df$peak_start+ann.df$summit,end=ann.df$peak_start+ann.df$summit,strand=ann.df$strand)},
	gr={ann=as.data.frame(ann.df)}
	)
	
	#ann: seqnames, start, end, strand
	Centres=ifelse(ann$strand==1|ann$strand=="+",ann$start,ann$end)
	if (anntype!="peak" & anntype!="annpeak") 
	{	ann$start=Centres+left
		ann$end=Centres+right
	}
	
	
	chr=as(ann$seqnames,"vector")
	start=mapply(max,ann$start,1)
	
	message("calculating")
	occupancy=c()
	for (i in 1:nrow(ann))
	{	
		cvg=cover.rpm[[chr[i]]]
		end=min(ann$end[i],length(cvg))
		if (start[i]<=end)
		{ 
			vtmp=cvg[start[i]:end]
			occupancy[i]=mean(vtmp)
			#occupancy[i]=sum(vtmp)
		}
	}
	
	tmp=ann.df
	#tmp$tss=Centres
	#tmp$up=ann$start
	#tmp$down=ann$end
	tmp$occupancy=occupancy

	tmp
}


####################### traveling_ratio ###################################	
traveling_ratio=function(ann.df,cover.rpm,left= -30,right=300)
{
	message("calculating promoter")
	ann=ann.df
	Centres=ifelse(ann$strand==1|ann$strand=="+",ann$start,ann$end)
	ann$start=Centres
	ann$end=Centres
	promoter=occupancy("tss",ann,cover.rpm,left,right)
	
	message("calculating body")
	ann=ann.df
	ann$start=ifelse(ann$strand==1|ann$strand=="+",ann$start+right,ann$start)
	ann$end=ifelse(ann$strand==1|ann$strand=="+",ann$end,ann$end-right)
	body=occupancy("peak",ann,cover.rpm,left,right)
	
	ann.df$promoter=promoter$occupancy
	ann.df$body=body$occupancy
	ann.df$traveling_ratio= promoter$occupancy / body$occupancy
	ann.df$traveling_ratio[is.na(ann.df$traveling_ratio)]=1
	
	ann.df

}