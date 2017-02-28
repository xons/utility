# customized IO
library("GenomicRanges")
library("GenomicAlignments")
library(parallel)
library(Rsamtools)

####################################################################################
####################################### Read in ####################################
df2gr=function(da,type='tss')
{
	library(GenomicRanges)
	switch(type,
		#annpeak output or as.data.frame(tr.gr)
	tss={gr=GRanges(da$seqnames,IRanges(da$start,da$end),strand=da$strand)},
		#ensembl data
	ensembl.tr={gr=GRanges(da$chromosome_name,IRanges(da$transcript_start,da$transcript_end),strand=da$strand)},
	ensembl.gene={gr=GRanges(seqnames=da$chromosome_name,IRanges(ann.df$start_position,ann.df$end_position),strand=da$strand)},
		#peak
	peak={gr=GRanges(da$chr,IRanges(da$start,da$end),strand=rep("*",nrow(da)),peak_index=rownames(da))},
	annpeak={gr=GRanges(da$seqnames,IRanges(da$peak_start,da$peak_end),strand=rep("*",nrow(da)),peak_index=da$peak_index)},
		#peak summit
	summit={da=GRanges(da$chr,IRanges(da$start+da$summit,da$start+da$summit),strand=rep(1,nrow(da)))},
	annsummit={da=GRanges(da$seqnames,IRanges(da$peak_start+da$summit,da$peak_start+da$summit),strand=da$strand)},
		#peak center
	center={da=GRanges(da$seqnames,IRanges(floor((da$start+da$end)/2),floor((da$start+da$end)/2)),strand=da$strand)},
	anncenter={da=GRanges(da$seqnames,IRanges(floor((da$peak_start+da$peak_end)/2),floor((da$peak_start+da$peak_end)/2)),strand=da$strand)},
		#others
	other={gr=GRanges(da[,1],IRanges(da[,2],da[,3]),strand=rep("*",nrow(da)))}	
	)
	#if ("summit" %in% names(da)) gr$summit=da$summit
	#seqlengths(gr)=seql
	gr
}

## read a chunck of bamfile
# bf=BamFile(fpath,yieldSize=100000)
# aln=readGAlignmentPairs(bf, param=param)
if (FALSE)
{#to save memory
bf=BamFile(fpath,yieldSize=100000)
open(bf,"r")
repeat{
 chunk=readGAlignmentPairs(bf,param=param)
 message(length(chunk))
 if (length(chunk) == 0L)
 break
	cvg.chunk = coverage(chunk)
	if (is.null(cvg)) {
	cvg = cvg.chunk
	} else {
	#library(pasillaBamSubset)
  	cvg = cvg + cvg.chunk
	}
}
close(bf)
}

readbam = function(fpath, paired=FALSE, unique.map=FALSE)
{
message("loading BAM file ...") 	
	if ( !paired & !unique.map ) 
	{
		flag=scanBamFlag(isUnmappedQuery = FALSE, isSecondaryAlignment = FALSE) ## get primary reads in muti-alignment(= -g 1 -x 1)
		#isDuplicate = FALSE ## remove duplicated reads
		param <- ScanBamParam(flag=flag, what="mapq",tag="NH")
		aln=readGAlignments(fpath, param=param)	
		message(length(aln))
	} 
	if ( !paired & unique.map )
	{
		flag=scanBamFlag(isUnmappedQuery = FALSE) 
		param <- ScanBamParam(flag=flag, what="mapq",tag="NH")
		aln=readGAlignments(fpath, param=param)
		message(paste0("mapped reads ",length(aln)))
		x=mcols(aln)
		aln=aln[x$NH==1] # get uniquely mapped reads in muti-alignment, if use -g 1 -x 1, this will not work, because all the reads are marked as NH=1
        message(paste0("uniquely mapped reads ",length(aln)))
	}
	if ( paired & !unique.map ) 
	{
		flag=scanBamFlag(isUnmappedQuery = FALSE,  isProperPair= TRUE, isSecondaryAlignment = FALSE) ## get 1 primary reads in muti-alignment(= -g 1 -x 1)
		param <- ScanBamParam(flag=flag, what="mapq",tag="NH")  ##??????  need isSecondaryAlignment = FALSE and tag="NH" ??????????
		aln=readGAlignmentPairs(fpath, param=param)#!! cannot get mapq from readGAlignmentPairs
		#x=mcols(aln)
		#aln=aln[x$mapq > 0] ## ?????for bwa alignment >10
		message(length(aln))
	}

#seqlevels(aln,force=TRUE) <- seqlevels(aln)[nchar(seqlevels(aln)) <10]

aln
}



nakename=function(fpath)
{
if ("/" %in% unlist(strsplit(as.vector(fpath[1]),""))) x=basename(fpath) else x=fpath
tools::file_path_sans_ext(x)
}

###### read in peaks and summits from MACS results 
read.macs=function(fname,fdr=5,p=50)
{
fdir=paste("MACS/",fname,"_peaks.xls",sep="")
peaks=read.table(fdir,comment="#",header=T,sep="\t")
peaks=subset(peaks,chr!="chrM")
if ("FDR..." %in% names(peaks))
subpeaks=subset(peaks, FDR... <= fdr & X.10.log10.pvalue.> p) else
subpeaks=subset(peaks, X.10.log10.pvalue.> p)
#subpeaks=subset(peaks, X.10.log10.pvalue.>30 & tags > 10 )
rownames(subpeaks)=1:nrow(subpeaks)
subpeaks

}

read.macs2=function(fname,p=5,reads=20)
{
#fdir=paste0("./",fname,"_peaks.xls")
peaks=read.delim(fname,comment="#")
peaks$chr=as.vector(peaks$chr)
peaks=subset(peaks,chr!="chrM")
peaks=peaks[!(1:nrow(peaks) %in% grep("chrUn",peaks$chr)),]
subpeaks=subset(peaks, X.log10.pvalue.> p & pileup > reads)
rownames(subpeaks)=1:nrow(subpeaks)
subpeaks

}

read.sicer=function(fdir,fdr=1e-10)
{
 peaks=read.table(fdir,comment="#",header=F,sep="\t")
names(peaks)=c("chr","start","end","counts","control","p_value","fold_change","FDR")
peaks=subset(peaks,chr!="chrM")
subpeaks=subset(peaks,FDR < fdr)
rownames(subpeaks)=1:nrow(subpeaks)
subpeaks
}

###### read files into RleLists
bw2coverage= function (fdir)
{
library(rtracklayer)
message('read bigwig file...')
x=import.bw(fdir)
#x=import.bw(fdir,selection = BigWigSelection(peak.gr))
coverage(x,weight=x$score)
}

wig2coverage=function(fdir)
{ 
library(rtracklayer)
df2Rle=function(df,step) 
{ 
	v0=df$score
	va=array(v0,c(length(v0),step))
	n=df$start[1]-step/2-1
	v=rep(0,n)
	v=c(v,va)
	Rle(v)
}

x=import(fdir,format="wig")
x.df=as.data.frame(x)
step=x.df$start[2]-x.df$start[1]
x.ls=split(x.df,x.df$space)
x.rle=lapply(x.ls,df2Rle,step)

}

## reform cuffdiff results
reform.exp=function(expkd)
{
name=deparse(substitute(expkd))
x=expkd[,c(2,3,7:10,12:14)]
names(x)[4:9]=c("wt",name,"log2fc","p","q","sig")
names(x)[c(3:4,6:9)]=paste(names(x)[c(3:4,6:9)],name,sep=".")
x
}


####################################################################################
##################################### Write out ####################################

##### for view in treeview, convert dataframe to cdt format
df2cdt=function(dataframe, fname)
{
dataframe=as.data.frame(dataframe)
NAME=rownames(dataframe)
GWEIGHT=rep(1,nrow(dataframe))
tmp=cbind(NAME,GWEIGHT,dataframe)
EWEIGHT=c(NA,NA,rep(1,ncol(dataframe)))
tmp=rbind(EWEIGHT,tmp)
rownames(tmp)[1]="EWEIGHT"
fdir=paste("cluster/",fname,".cdt",sep="")
#fdir=paste0(fname,".cdt")
write.table(tmp,file=fdir,sep="\t",col.names=NA)
}

######
write.na=function(x,fdir,add=F)	#write vannila, no colname, no rowname
{
	write.table(x,fdir,sep="\t",col.names=F,row.names=F,quote=F,append=add)
}
write.nsv=function(x,fdir,add=F)	#write vannila
{
	write.table(x,fdir,sep=",",col.names=F,row.names=F,quote=F,append=add)
}
write.cn=function(x,fdir,add=F)	#write with column names
{
	write.table(x,fdir,sep="\t",col.names=T,row.names=F,quote=F,append=add)
}

write.txt=function(x,fname, rowname.header=" ") #write with colname and rowname, give a header for rowname column
{
 #fdir=paste0("cluster/",fname,".txt")
 fdir=paste0(fname,".txt")
 write.table(rowname.header,fdir,eol="\t",col.names=F,row.names=F,quote=F)
 write.table(x,fdir,sep="\t",col.names=T,row.names=T,quote=F,append=T)
}

write.gr2bed=function(gr,fname)
{
	df=as.data.frame(gr)
	write.tsv(df[,1:3],fname)
}
