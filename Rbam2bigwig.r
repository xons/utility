# Rscript ~/cmd/Rbam2bigwig.r 
# R script
# convert bam file to bigwig file
# cd to the output folder, if need the track line output

# fpath=($(ls *.bam))

library(optparse)

option_list = list(
	make_option("--seq", type="character", default=NULL, 
              help="Sequencing type, can be chip, rna or proseq", metavar="character"),
	make_option("--input", type="character", default=NULL, 
              help="Space separated Bam file paths, if more than one, use quote", metavar="character"),
	make_option("--stranded", type="logical", default=TRUE, 
              help="Stranded library or not, for rna Sequencing type only. [default= %default]", metavar="logical"),
	make_option("--paired", type="logical", default=FALSE, 
              help="Paired end or not. [default= %default]", metavar="logical"),
	make_option("--frag_min", type="integer", default=1, 
              help="Minimum fragment size. [default= %default]", metavar="integer"),
	make_option("--frag_max", type="integer", default=1e6, 
              help="Maximum fragment size. [default= %default]", metavar="integer")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser)
print(opt)
opt$input=unlist(strsplit(opt$input,split=" "))

## args <- commandArgs(TRUE)
## fnames.vector=args[c(-1)]

library(parallel)
library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(Rsamtools)
source("~/cmd/customIO.r") ## read bam


gr2bw = function(gr,fbw,n=NA,desc=NULL,color="0,0,255")
{ ## convert coverage GRanges to bigwig
message("calculating coverage ...")
if (is.na(n)) n=length(gr)
rpm=coverage(gr)/n*1000000
#tmp=aggregate(rpm, FUN = mean, start = seq(1,length(rpm)-25,by=25), width=25) # bin by window 25
#rpm=tmp

message("exporting bigWig file ...")
export.bw(rpm,fbw)


#### write UCSC track upload URL
fname=basename(fbw)
#desc=ifelse(paired,' FPM',' RPM')
#desc=ifelse(paired,' FPM',' RPM ext150')

if (n<0) color='255,0,0'

subpath=unlist(strsplit(getwd(),"/"))
trackpath=paste0("http://myserver/",subpath[8],"/",subpath[10])
trackfile="track.txt"

trackline=paste0('track type=bigWig name=',fname,' description=',fname, desc, ' graphType=bar visibility=full color=',color,' itemRgb=on autoScale=on bigDataUrl=',trackpath,"/",fname,'.bw')
write.table(trackline,trackfile,col.names=F,row.names=F,quote=F,sep="\n",append=T)
}

## for CHIPseq
Rbam2bigwig = function (fpath,paired=FALSE, frag_min=1, frag_max=1e6)
{
aln=readbam(fpath,paired=paired)

gr=granges(aln)

message("converting to GRanges ...")
if (frag_min!=1 | frag_max!=1e6)
{
##size select
x=width(gr) 
#z=quantile(x[x>80],c(.025,0.975))#116   218
#maln.gr=maln.gr[x>z[1] & x<z[2]]
gr=gr[x>frag_min & x<frag_max]
}

if (!paired)
{
tmp=resize(gr,150) # extend the reads to 150 bp towards middle
gr=trim(tmp)
}

desc=ifelse(paired,' FPM',' RPM ext150')
fbw=sub("bam|sorted.bam","bw",fpath)
gr2bw(gr,fbw,desc=desc)
}

## for RNAseq
Rbam2bigwig_RNA = function (fpath,stranded=FALSE, paired=FALSE)
{
aln=readbam(fpath,paired=paired)

message("converting to GRanges ...")
grl=grglist(aln)
gr=unlist(grl)

desc=ifelse(paired,' FPM',' RPM')
fbw=sub(".bam|.sorted.bam|.accepted_hits.bam",".bw",fpath)
gr2bw(gr,fbw,desc=desc)

if (stranded) 
{
n=length(gr)
gr1=gr[strand(gr)=="-"]
gr2=gr[strand(gr)=="+"]

fbw1=sub(".bam|.sorted.bam|.accepted_hits.bam","_positiveStrand.bw",fpath)
fbw2=sub(".bam|.sorted.bam|.accepted_hits.bam","_negativeStrand.bw",fpath)
gr2bw(gr1,fbw1,n,desc=desc)
gr2bw(gr2,fbw2,-n,desc=desc)
}

}

Rbam2bigwig_Proseq = function (fpath,paired=FALSE)
{
aln=readbam(fpath,paired=paired)

message("converting to GRanges ...")
maln.gr=granges(aln)

maln.end.gr=resize(maln.gr,1,fix='end') ## use the last bp only

n=length(maln.gr)
maln1.gr=maln.end.gr[strand(maln.end.gr)=="+"] ## no need to switch strand if did reverse complement for fastq
maln2.gr=maln.end.gr[strand(maln.end.gr)=="-"]

fbw1=sub(".bam|.sorted.bam","_positiveStrand.bw",fpath)
fbw2=sub(".bam|.sorted.bam","_negativeStrand.bw",fpath)
gr2bw(maln1.gr,fbw1,n)
gr2bw(maln2.gr,fbw2,-n)

}


switch( opt$seq,	
chip={ mclapply(opt$input, Rbam2bigwig, paired=opt$paired, frag_min=opt$frag_min, frag_max=opt$frag_max, mc.cores = 6) },
rna={ mclapply(opt$input, Rbam2bigwig_RNA, stranded=opt$stranded, paired=opt$paired, mc.cores = 6) },
proseq={ mclapply(opt$input, Rbam2bigwig_Proseq, paired=opt$paired, mc.cores = 6) }
)




########## old #############
## for stranded RNAseq
Rbam2bigwig_RNA_stranded = function (fdir,paired=FALSE,tophat=TRUE)
{
message("loading BAM file ...")
aln=readbam(fpath,paired=paired,tophat=tophat)

message("converting to GRanges ...")
#n=length(aln)
#aln1=aln[strand(aln)=="-",]
#aln2=aln[strand(aln)=="+",]

#aln1.grl=grglist(aln1)
#maln1.gr=unlist(aln1.grl)
#aln2.grl=grglist(aln2)
#maln2.gr=unlist(aln2.grl)

aln.grl=grglist(aln)
maln.gr=unlist(aln.grl)
n=length(maln.gr)
maln1.gr=maln.gr[strand(maln.gr)=="-"]
maln2.gr=maln.gr[strand(maln.gr)=="+"]
rm(maln.gr)

fbw1=paste0(fdir,"_positiveStrand.bw")
fbw2=paste0(fdir,"_negativeStrand.bw")
gr2bw(maln1.gr,fbw1,n)
gr2bw(maln2.gr,fbw2,n,reverse=TRUE)

}


## for paired-end chip-seq/ DNAseq
Rbam2bigwig_pair = function (fdir)
{
message("loading BAM file ...")
flag=scanBamFlag(isUnmappedQuery = FALSE,  isProperPair= TRUE) 
#flag=scanBamFlag(isUnmappedQuery = FALSE, isDuplicate = FALSE) ## remove duplicated reads
param <- ScanBamParam(flag=flag, what="mapq")
aln=readGAlignmentPairs(fpath, param=param)
#seqlevels(aln,force=TRUE) <- seqlevels(aln)[nchar(seqlevels(aln)) <10]
#x=mcols(aln)
#aln=aln[x$mapq > 0]   

maln.gr=granges(aln)
rm(aln) #save memory

message("calculating coverage ...")
maln.rpm=coverage(maln.gr)/length(maln.gr)*1000000

message("exporting bigWig file ...")
fbw=sub("bam|sorted.bam","bw",fpath)
export.bw(maln.rpm, fbw)

#write UCSC track upload URL
fname=basename(fbw)
trackline=paste0('track type=bigWig name=',fname,' description="',fname,' FPM" graphType=bar visibility=full color="0,0,255" itemRgb=on autoScale=on bigDataUrl=',trackpath,"/",fname)
write.table(trackline,trackfile,col.names=F,row.names=F,quote=F,sep="\n",append=T)
#write.table(t(c("MNase",subpath[c(6:8,10)],fname,tmp)),tracktable,col.names=F,row.names=F,quote=F,sep="\t",append=T)

}


## for 4C seq
Rbam2bigwig_4c = function (fname)
{
#folder="kca/GSE37275"
message("loading BAM file ...")
fpath=paste(fdir,"/",fname,".sorted.bam",sep="")
aln=readGAlignments(fpath)
#n=length(aln)
#seqlevels(aln,force=TRUE) <- seqlevels(aln)[nchar(seqlevels(aln)) <10]
message("converting to GRanges ...")
maln.gr=granges(aln)

libraryName= "hind3.dpn2.mm9.51bp.csv"
lib=read.csv(paste0("4c/",libraryName),sep="\t")

#lib.gr=peak.gr(lib)
#maln.sub.gr=subsetByOverlaps( maln.gr,lib.gr, maxgap = 0L, minoverlap = 1L, type ="any", ignore.strand = T)
 lib.leftend=lib[,1:3]
 lib.leftend[,3]=lib.leftend[,2]+lib$leftFragEndLength
 lib.rightend=lib[,1:3]
 lib.rightend[,2]=lib.rightend[,3]-lib$rightFragEndLength
 lib.end.gr=peak.gr(rbind(lib.leftend,lib.rightend))
 maln.sub.gr=subsetByOverlaps( maln.gr,lib.end.gr, maxgap = 0L, minoverlap = 1L, type ="any", ignore.strand = T)

#tmp=resize(maln.gr,150) # extend the reads 150 bp towards middle
#maln.gr.extend=tmp

message("calculating coverage ...")
maln.rpm=coverage(maln.sub.gr)/length(maln.sub.gr)*1000000

message("exporting bigWig file ...")
fbw=paste(fdir,"/",fname,".bw",sep="")
export.bw(maln.rpm, fbw)

#write UCSC track upload URL
tmp=paste('track type=bigWig name=',fname,' description="',fname,' rpm" graphType=bar visibility=full color="0,0,255" itemRgb=on autoScale=on bigDataUrl=http://myserver',trackpath,"/",fname,'.bw',sep="")
write.table(tmp,trackline,col.names=F,row.names=F,quote=F,sep="\n",append=T)
write.table(t(c("4C",subpath[c(6:8,10)],fname,tmp)),tracktable,col.names=F,row.names=F,quote=F,sep="\t",append=T)
}

