########### for gene expression
library("GenomicRanges")
library(GenomicAlignments)
library(parallel)
library(Rsamtools)
source("customIO.r") ## read bam

count.gene=function(fnamels,genome,output,stranding=TRUE,paired=FALSE,unique.map=FALSE,antisense=FALSE)
{
    genome=paste("gnModel",genome,sep=".")
    load(genome)
    countls=mclapply(fnamels,count.g,gnModel,stranding,paired,unique.map,antisense,mc.cores=6,mc.set.seed = FALSE)
    #names(countls)=tools::file_path_sans_ext(basename(fnamels))
	fnamels=sub(".bam|.sorted.bam|.accepted_hits.bam", "",fnamels)
	names(countls)=basename(fnamels)
	#save(countls,file="tmp")
	
    count.table= do.call(cbind,lapply(countls,"[[",1))
    size=unlist(lapply(countls,"[[",2))  
	count.norm=t(t(count.table)/size)*1e6
	#colnames(count.norm)=paste("norm",colnames(count.table),sep=".")
	if (exists('exon.length')==TRUE) count.fpkm=count.norm/exon.length*1000
	#colnames(count.fpkm)=paste("fpkm",colnames(count.table),sep=".")
	#count.log2fpkm=log2(count.fpkm)
	#colnames(count.log2fpkm)=paste("log2fpkm",colnames(count.table),sep=".")
	
    save(countls,count.table,size,count.fpkm,file=paste0(output,".exp.countls"))
	
	png(paste0(output,".all.correlation.png"),width=480*max(c(1,ncol(count.table)/50)),height=480*max(c(1,ncol(count.table)/50)))
		count.log2fpkm=log2(count.fpkm)
		corplot(count.log2fpkm)
	dev.off()
    }

count.g=function(fpath,gnModel,stranded=TRUE,paired=FALSE,unique.map=FALSE,antisense=FALSE)
{	
	aln=readbam(fpath,paired=paired, unique.map=unique.map)
	
	n=length(aln)
	print(n)
	if ( !stranded ) strand(aln) <- "*" # for strand-blind sample prep protocol
	if ( stranded & !antisense )
    {# for strand specific, need to flip the strand
	print("stranded")
	tmp1=which(strand(aln)=="+")
	tmp2=which(strand(aln)=="-")
	strand(aln)[tmp1]=Rle("-",length(tmp1))
	strand(aln)[tmp2]=Rle("+",length(tmp2))
	}
		#counts=summarizeOverlaps(gnModel, aln, mode="IntersectionNotEmpty")
#	hits <- countOverlaps(aln, gnModel) # find multi-alignment reads
#    counts <- countOverlaps(gnModel, aln[hits==1]) # count without reads cover more than one gene
	counts <- countOverlaps(gnModel, aln)
    names(counts) <- names(gnModel)
	counts[counts==0]=1
	return(list(counts,n))
}

corplot=function(m,main="Correlation of log2 FPKM",ann=NA,anncolor=NA)
{## plot correlation matrix
library(lattice)
rgb.palette <- colorRampPalette(c("blue", "yellow"), space = "rgb")
x=cor(m)

# method 1
#hc=hclust(dist(x))
#levelplot(x[hc$order,hc$order], main="correlation matrix", xlab="", ylab="", col.regions=rgb.palette(120), cuts=100, at=seq(0,1,0.01)) 
#levelplot(x[hc$order,hc$order], main="correlation matrix", xlab="", ylab="", col.regions=rgb.palette(120))

# method 2
library(pheatmap)
pheatmap(x,col=rgb.palette(120),main=paste0(main," N=", nrow(m)), annotation_row = ann, annotation_colors =anncolor)

# method 3
#heatmap(cor(m),col=rgb.palette(120))

# method 4
#library(gplots)
#heatmap.2(cor(m),col=rgb.palette(120),trace="none",density.info="none",keysize=1)
}
	
############### for peak diff test
count.peak=function(fnamels,peak_gr,output,paired=FALSE,stranded=FALSE,antisense=FALSE,extend=TRUE)
{
	peak_gr
	countls=mclapply(fnamels,count.p,peak_gr,paired=paired, stranded=stranded, antisense=antisense, extend=extend, mc.cores=6,mc.set.seed = FALSE)
	#names(countls)=tools::file_path_sans_ext(basename(fnamels))
	fnamels=sub(".bam|.sorted.bam|.accepted_hits.bam", "",fnamels)
	names(countls)=basename(fnamels)
    count.table= do.call(cbind,lapply(countls,"[[",1))
    size=unlist(lapply(countls,"[[",2))   
	count.norm=t(t(count.table)/size)*1e6
	colnames(count.norm)=paste("norm",colnames(count.norm),sep=".")
	count.fpkm=count.norm/width(peak_gr)*1000
	colnames(count.fpkm)=paste("fpkm",colnames(count.table),sep=".")
	save(countls,count.table,size,count.norm, count.fpkm, file=paste0(output,".chip.countls"))
	
	count.log2fpkm=log2(count.fpkm)
	colnames(count.log2fpkm)=paste("log2fpkm",colnames(count.table),sep=".")
	# save(countls,count.table,size,count.norm,count.fpkm,count.log2fpkm,file=paste0(output,".exp.countls"))
	if (length(fnamels)>1)
	{
	png(paste0(output,".all.correlation.png"))
		colnames(count.log2fpkm)=colnames(count.table)
		corplot(count.log2fpkm)
	dev.off()
	}
	
}

count.p=function(fpath,peak_gr,paired=FALSE,stranded=FALSE,antisense=FALSE,extend=TRUE)
{
	aln=readbam(fpath,paired=paired)
	n=length(aln)
	message(n)

	if (stranded==F) 	strand(aln) <- "*" # for strand-blind sample prep protocol
	if (stranded==T & antisense==F)
    {# for strand specific, need to flip the strand
	print("stranded")
	tmp1=which(strand(aln)=="+")
	tmp2=which(strand(aln)=="-")
	strand(aln)[tmp1]=Rle("-",length(tmp1))
	strand(aln)[tmp2]=Rle("+",length(tmp2))
	}
	
	message("coverting to GRanges")
	maln.gr=granges(aln)
	
	if (paired==FALSE & extend==TRUE) 
	{	
	tmp=resize(maln.gr,150) # extend the reads 150 bp towards middle
	maln.gr=trim(tmp)
	}
	message("counting")
		#hits <- countOverlaps(aln, peak_gr)
		#counts <- countOverlaps(peak_gr, aln[hits==1])	#discards a read hitting more than one peak
    	#countswt <- countOverlaps(peak_gr, alnwt)
	counts <- countOverlaps(peak_gr, maln.gr)
	counts[counts==0]=1
	names(counts)=names(peak_gr)
	return(list(counts,n))
}

