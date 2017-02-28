#source("/n/projects/xig/cmd/metagene.r")
# mclapply(rpmls,metagene,ann.df=annmll_19_loss,mc.cores=6,mc.set.seed = FALSE)

metagene = function(  query, ann, anntype="tss",winstart=c(300,300),winend=c(300,300), bin_width=10, bodylength=1000 )
{
# query can be file path of bigwig or coverage of RleList
# ann can be data frame or GRanges
# win is the bp on each side
	library(GenomicRanges)
	if (class(query)=="character") cover.rpm=bw2coverage(query) else cover.rpm = query
	
	if (class(ann)=='data.frame') ann.gr=df2gr(ann,type=anntype) else ann.gr=ann
	seql=unlist(lapply(cover.rpm,length))
	seqlengths(ann.gr)=seql[seqlevels(ann.gr)]
	
#	tss.gr=resize(ann.gr,1,fix="start")
#	left=resize(tss.gr,winstart[1]+1,fix="end");left=trim(left) 
#	right=resize(tss.gr,winstart[2]+1,fix="start");right=trim(right)
#	start.gr=punion(left,right)	
	start.gr=promoters(ann.gr,winstart[1],winstart[2])
	start.gr=trim(start.gr)
	
	tes.gr=resize(ann.gr,1,fix="end")
#	left=resize(tes.gr,winend[1]+1,fix="end");left=trim(left)
#	right=resize(tes.gr,winend[2]+1,fix="start");right=trim(right)
#	end.gr=punion(left,right)
	end.gr=promoters(tes.gr,winend[1],winend[2])
	end.gr=trim(end.gr)
	
	body1.gr=psetdiff(ann.gr,start.gr)
	body.gr=psetdiff(body1.gr,end.gr)

	nrows=length(ann.gr)
	v = array( 0, c(nrows,floor((sum(winstart,winend)+bodylength)/10)))
	
	message("calculating")	
	s.ls=cover.rpm[start.gr]
	e.ls=cover.rpm[end.gr]
	b.ls=cover.rpm[body.gr]
	#b.norm=lapply(b,function(x){Rle(approx(x,n=bodylength)$y)})
	#b.norm=as(b.norm,'RleList')
	
	for( i in 1:nrows ) 
	{	
#		s=cover.rpm[start.gr[i]][[1]]
#		b=(approx(cover.rpm[body.gr[i]][[1]],n=bodylength))$y
#		e=cover.rpm[end.gr[i]][[1]]

		s=s.ls[[i]]
		e=e.ls[[i]]
		b=(approx(b.ls[[i]],n=bodylength))$y
		
		strand=as.vector(strand(start.gr[i]))
		if( strand== -1| strand=="-" ) vtmp=rev(c(as.vector(e),b,as.vector(s))) else
		vtmp=c(as.vector(s),b,as.vector(e))
		
		if (width(start.gr[i])< sum(winstart)) vtmp=c(rep(0,sum(winstart)-width(start.gr[i])),vtmp)
		if (width(end.gr[i])< sum(winend)) vtmp=c(vtmp,rep(0,sum(winend)-width(end.gr[i])))
		
		v[i,]=bin_mean(vtmp,bin_width)	
 	}
	v
}

bin_mean=function(v,bin)
{	
	colMeans(array(v,c(bin,length(v)/bin)))
}

ggplot.meta=function(mxls,fname,gene.group,winstart=c(300,300),winend=c(300,300),bin=10,bodylength=1000,ylim=NA)
{
## gene.group is a factor of group with gene id as names, same order as mx 
  library(reshape2)
  library(ggplot2)
  nls=lapply(mxls,function(x){rownames(x)=names(gene.group);x})
  da=melt(nls)
  colnames(da)=c('id','bin','coverage','samples')
  da$group=gene.group[match(da$id,names(gene.group))]  

    g=ggplot(da,aes(bin,coverage,color=samples))+
          geom_vline(xintercept = c(winstart[1],bodylength+sum(winstart)+winend[1])/bin,color='gray60',linetype='solid',size=0.5)+
          geom_vline(xintercept = c(sum(winstart),bodylength+sum(winstart))/bin,color='gray60',linetype='longdash',size=0.5)+
          geom_smooth()+ ## Gray region is the 95% confidence level interval for predictions from a generalized additive model. ( +/- 1.96*standard deviation )
          scale_x_continuous(breaks=c(0,winstart[1],sum(winstart),sum(winstart)+bodylength/2,bodylength+sum(winstart),bodylength+sum(winstart)+winend[1],bodylength+sum(winstart)+sum(winend))/bin
				,labels=c(paste0('-',winstart[1]),'TSS',paste0('+',winstart[2]),'scaled body',paste0('-',winend[1]),'TES',paste0('+',winend[2]))) +        
          facet_grid(. ~ group)+
          ggtitle(fname)+
          theme(text=element_text(size=15),plot.title = element_text(hjust = 0.5))+
          xlab('')+
          ylab('coverage (RPM)')
	if (!is.na(ylim)) g=g+coord_cartesian(ylim=ylim)
    ggsave(paste0('metaplot.',fname,'.png'),g,width=7*length(levels(da$group)))  
  }

#  f(mxls,'test')

  
plot.meta=function(meta.out,winstart=c(300,300),winend=c(300,300),bin=10,ymax=4,col="black",add=F,bodylength=1000)
{
b=bodylength/bin
if(add==T) lines(colSums(meta.out,na.rm=T)/nrow(meta.out),lwd=3,col=col) else
	{
plot(colSums(meta.out,na.rm=T)/nrow(meta.out),type="l",lwd=3,col=col,ylim=c(0,ymax),axes=F)
abline(v=winstart[1]/bin,col="gray",lwd=3) #TSS
abline(v=(sum(winstart)+winend[1])/bin+b,col="gray",lwd=3) #TES
abline(v=sum(winstart)/bin,col="gray")
abline(v=sum(winstart)/bin+b,col="gray")
axis(2,at=seq(0,ymax,length.out=5),labels=seq(0,ymax,length.out=5))
axis(1,at=c(winstart[1]/bin,0,sum(winstart)/bin,(sum(winstart)+winend[1])/bin+b,sum(winstart)/bin+b,sum(winstart,winend)/bin+b),labels=c("TSS",-winstart[1],winstart[2],"TES",-winend[1],winend[2]))
	}
 #savePlot("meta.png")
 }

 plot.meta.median=function(meta.out,winstart=c(300,300),winend=c(300,300),bin=10,ymax=4,col="black",add=F,bodylength=1000)
{
b=bodylength/bin
cols=col2rgb(col)
tmp=apply(meta.out,2,fun)

if(add==F) 
   {
plot(tmp["median",],type="l",lwd=3,col=col,ylim=c(0,ymax),axes=F)
abline(v=winstart[1]/bin,col="gray",lwd=3) #TSS
abline(v=(sum(winstart)+winend[1])/bin+b,col="gray",lwd=3) #TES
abline(v=sum(winstart)/bin,col="gray")
abline(v=sum(winstart)/bin+b,col="gray")
axis(2,at=seq(0,ymax,length.out=5),labels=seq(0,ymax,length.out=5))
axis(1,at=c(winstart[1]/bin,0,sum(winstart)/bin,(sum(winstart)+winend[1])/bin+b,sum(winstart)/bin+b,sum(winstart,winend)/bin+b),labels=c("TSS",-winstart[1],winstart[2],"TES",-winend[1],winend[2]))
    } 
polygon(c(1:dim(tmp)[2],rev(1:dim(tmp)[2])),c(tmp["confu",],rev(tmp["confl",])),col=rgb(red=cols[1], green=cols[2], blue=cols[3],alpha=32,maxColorValue=255),border=NA)
#plot(tmp["confu",],type="l",lty=2)
lines(tmp["median",],type="l",lwd=3,col=col)
#lines(tmp["confl",],type="l",lty=2)

 }

 
 #f=function(v) {quantile(v,c(.25, 0.5,0.75))}
fun=function(v)
{
	tmp=boxplot.stats(v)
     y=c(tmp$stats[2:4],tmp$conf)
     names(y)=c("25","median","75","confl","confu")
     y
}



