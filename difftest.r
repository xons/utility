edger.test=function(input,column,group,genome,output,norep=F,batch=c(),condition,fdr=0.05, logfc=0,gomap="org.Dm.eg.db",size.new=NA, norm='TMM')
{
 #input="sox2mutant"
 #column=c(7:9,1:6)
 #group=c(1,1,1,2,2,2,3,3,3)
 #genome="mm10"
 #output="sox2.batch"
 #batch=c(1,2,3,1,2,3,1,2,3)
 #norep=F
 #condition=c("s2A","s2E") # name of the test samples

annpath='annotation/'
####### test #########
message("test") 
    library(edgeR)  
    load(paste0(input,".exp.countls"))      
 
    counts=count.table[,column]
	#log2fpkm=count.log2fpkm[,column]
	cutoff=1 
	index=which(rowSums(counts >cutoff)>=2)
	#index=which(rowSums(cpm(counts)>cutoff)>=2)
	#index=which(rowSums(count.fpkm[,column] > cutoff) >= 2)
	counts=counts[index,]
	
	group=as.factor(group)
	dge <- DGEList(counts, group=group)
    if (!is.na(size.new)) dge$samples$lib.size=size.new[column]
	print(norm)
    dge <- calcNormFactors(dge,method=norm)
	
	if (!is.null(batch))
	{   design <- model.matrix( ~ group + batch, data=dge$samples ) ## removeBatchEffect(x,batch) # x must be log expression value
	} else
		design <- model.matrix( ~ group, data=dge$samples )
	
	load(paste0(annpath,"gnModel.",genome))
	fpkms=rpkm(dge,gene.length=exon.length[rownames(dge)])
	
	png(paste0(output,".exp.correlation.png"))
	#log2fpkm=log2fpkm[index,]	
	#colnames(log2fpkm)=colnames(counts)
	log2fpkm=log2(fpkms)
	corplot(log2fpkm)
	if (!is.null(batch)) 
	{log2fpkm.b=removeBatchEffect(log2fpkm,batch=batch,design=design)
	corplot(log2fpkm.b)}
	dev.off()
	
	if (norep ==F) #plot MDS
	{
		png(paste0(output,".MDS.png"))
		plotMDS.DGEList(dge,col=as.numeric(group), method="bcv",main="MDS plot") #An MDS plots shows distances, in terms of biological coefficient of variation (BCV), between samples:
		dev.off()
	}

	if (norep==F) #fit
	{
		dge <- estimateGLMCommonDisp(dge, design) #When estimating tagwise dispersions, the empirical Bayes method is applied to squeeze tagwise dispersions towards common dispersion
													#or trended dispersions, whichever exists. If both exist, the default is to use the trended dispersions.
		dge <- estimateGLMTrendedDisp(dge, design)
		# dge <- estimateGLMTrendedDisp(dge, design, method="power")
		dge <- estimateGLMTagwiseDisp(dge, design)
		png(paste0(output,".BCV.png"))
		plotBCV(dge, main="BCV plot") #biological coeffcient of variation
		dev.off()
		fit <- glmFit(dge, design, prior.count=0)
	} else
		fit <- glmFit(dge, design, prior.count=0,dispersion=0.01 ) #0.4 for human data, 0.1 for data on genetically identical model organisms or 0.01 for technical replicates.
		
########### output ###########	
## output formatting
	n=length(condition)
		#lrTest <- glmLRT(fit, coef=2)
	lrTest=lapply(as.list(seq(2,n+1)),glmLRT,glmfit=fit)
		#tt=lrTest$table
	ttls=lapply(lrTest,"[[","table")
	names(ttls)=condition
	tmp=lapply(ttls,function(x){
								x$FDR=p.adjust(x$PValue, method="BH")
								return(x)
								})
	ttls=tmp
	tt=do.call(cbind,ttls)
		#tt$FDR=p.adjust(tt$PValue,method="BH")
		#tt=topTags(lrTest,Inf)
		#tt=lapply(lrTest,topTags,Inf)

	load(paste0(annpath,genome,".rda"))
	#genome=eval(parse(text=genome))
	#tmp=match(rownames(counts),genome$Ensembl.Gene.ID)
	tmp=match(rownames(counts),genome$ensembl_gene_id)
	#out.edger=data.frame(genome[tmp,],gene.id=rownames(counts),counts,count.fpkm[index,column],tt)
	out.edger=data.frame(genome[tmp,],gene.id=rownames(counts),fpkms,tt)
	# out.edger=out.edger[!is.na(out.edger$Chromosome.Name),]
	save(dge,fit,lrTest,ttls,out.edger,fpkms,file=paste0(output,".exp.edger")) 
	
#plot MA
message("plot MA")
	mapply(plotma_edger,ttls,names(ttls),fdr,logfc=logfc)	
	
#output tsv formatting
message("output table")
	tmp=lapply(ttls,function(x) {x=x[,c(2,1,4,5)]})
	ttselect=do.call(cbind,tmp)
	#out.tab=data.frame(genome[match(rownames(ttselect),genome$Ensembl.Gene.ID),c(1:13)],fpkms,ttselect)#reorder columns of tt
	out.tab=data.frame(genome[match(rownames(ttselect),genome$ensembl_gene_id),c(1:2,4:7,13,16)],fpkms,ttselect)#reorder columns of tt
	## xx=lapply(ttls,function(x) {x=which(x[4]<0.05)})
	## x=unique(unlist(xx))
	## out.tab=out.tab[x,] # select lines with either condition p<0.05
	## out.tab=out.tab[order(ttselect[x,3]),]# sort by the first pvalue
	out.tab[,-1:-8]=signif(out.tab[,-1:-8],3) # round the number
	write.cn(out.tab,paste0(output,".exp.edger.tsv"))	

## gene list
message("save gene list")
	mapply(get_genelist,ttls,names(ttls),MoreArgs=list(fdr=fdr,log2fc=logfc,genome=genome))
		
## go analysis with topGO
message("topgo")
	mapply(go.topgo,ttls,names(ttls),MoreArgs=list(score.cutoff=fdr,log2fc=0,mapping=gomap,genome=genome))

}

edger.contrast=function(input,column,group,genome,output,norep=F,batch=c(),contrast,fdr=0.05, logfc=0,gomap="org.Dm.eg.db")
{
 #input="sox2mutant"
 #column=c(7:9,1:6)
 #group=c(1,1,1,2,2,2,3,3,3)
 #genome="mm10"
 #output="sox2.batch"
 #batch=c(1,2,3,1,2,3,1,2,3)
 #norep=F
 #contrast=list(BcatKO.24h=c(-1,1,0,0,0,0,0),
 #             BcatKO.48h=c(0,0,-1,1,0,0,0),
 #             shBcat=c(0,0,0,0,-1,1,0),
 #             shE2F6=c(0,0,0,0,-1,0,1))

 annpath='annotation/'

 message("test") 
    library(edgeR)  
    load(paste0(input,".exp.countls"))      
 
    counts=count.table[,column]
	#log2fpkm=count.log2fpkm[,column]
	cutoff=1 #cutoff is 1 fpkm
	index=which(rowSums(counts >cutoff)>=2)
	#index=which(rowSums(cpm(counts)>cutoff)>=2)
	#index=which(rowSums(count.fpkm > cutoff) >= 2)
	counts=counts[index,]
		
	group=as.factor(group)
	dge <- DGEList(counts, group=group)
          #dge$samples$lib.size=size[column]
    dge <- calcNormFactors(dge)
		
	load(paste0(annpath,"gnModel.",genome))
	fpkms=rpkm(dge,gene.length=exon.length[rownames(dge)])
	
	png(paste0(output,".exp.correlation.png"))
	#log2fpkm=log2fpkm[index,]	
	#colnames(log2fpkm)=colnames(counts)
	log2fpkm=log2(fpkms)
	corplot(log2fpkm)
	dev.off()

	if (norep ==F) #plot MDS
	{
		png(paste0(output,".MDS.png"))
		plotMDS.DGEList(dge,col=as.numeric(group), method="bcv",main="MDS plot")#An MDS plots shows distances, in terms of biological coefficient of variation (BCV), between samples:
		dev.off()
	}

	if (!is.null(batch))
	{   design <- model.matrix( ~ 0+ group + batch, data=dge$samples )
	} else
		design <- model.matrix( ~ 0+group, data=dge$samples )

	if (norep==F) #fit
	{
		dge <- estimateGLMCommonDisp(dge, design)
		dge <- estimateGLMTrendedDisp(dge, design)
		dge <- estimateGLMTagwiseDisp(dge, design)
		png(paste0(output,".BCV.png"))
		plotBCV(dge,main="BCV plot") #biological coeffcient of variation
		dev.off()
		fit <- glmFit(dge, design, prior.count=0)
	} else
		fit <- glmFit(dge, design, prior.count=0,dispersion=0.01 )
#save(fit,file="tmp")
## output formatting
	f=function(contrast)
	{
		lrt = glmLRT(fit, contrast=contrast)
		x=lrt$table
		x$FDR=p.adjust(x$PValue, method="BH")
		x
	}	
	ttls=lapply(contrast,f)
	names(ttls)=names(contrast)
	tt=do.call(cbind,ttls)
		#tt=topTags(lrTest,Inf)
		#tt=lapply(lrTest,topTags,Inf)
	
	load(paste0(annpath,genome,".rda"))
	#genome=eval(parse(text=genome))
	#tmp=match(rownames(counts),genome$Ensembl.Gene.ID)
	tmp=match(rownames(counts),genome$ensembl_gene_id)
	out.edger=data.frame(genome[tmp,],gene.id=rownames(counts),fpkms,tt)
	save(dge,fit,ttls,out.edger,fpkms,file=paste0(output,".exp.edger")) 
	
	#plot MA
message("plot MA")
	mapply(plotma_edger,ttls,names(ttls),fdr,logfc=logfc)	
	
#output tsv formatting
message("output table")
	xx=lapply(ttls,function(x) {x=x[c(2,1,4,5)]})
	ttselect=do.call(cbind,xx)
	#out.tab=data.frame(genome[match(rownames(ttselect),genome$Ensembl.Gene.ID),c(1:13)],fpkms,ttselect)#reorder columns of tt
	out.tab=data.frame(genome[match(rownames(ttselect),genome$ensembl_gene_id),c(1:2,4:7,13,16)],fpkms,ttselect)#reorder columns of tt
	out.tab[,-1:-8]=signif(out.tab[,-1:-8],3) # round the number
	write.cn(out.tab,paste0(output,".exp.edger.tsv"))	

## gene list
message("save gene list")
	mapply(get_genelist,ttls,names(ttls),MoreArgs=list(fdr=fdr,log2fc=logfc,genome=genome))
		
## go analysis with topGO
message("topgo")
	mapply(go.topgo,ttls,names(ttls),MoreArgs=list(score.cutoff=fdr,log2fc=logfc,mapping=gomap,genome=genome))


}

get_genelist=function(tt,name,fdr,log2fc=0,genome)
{
	tt.up=subset(tt,FDR<fdr & logFC > log2fc)
	tt.down=subset(tt,FDR<fdr & logFC < -log2fc)
	
	#tmp=match(rownames(tt.up),genome$Ensembl.Gene.ID)
	tmp=match(rownames(tt.up),genome$ensembl_gene_id)
	gene.up=data.frame(genome[tmp,c(1:2,4:7,13,16)],tt.up)
	write.cn(gene.up,paste0(name,".genelist.up.tsv"))
	
	#tmp=match(rownames(tt.down),genome$Ensembl.Gene.ID)
	tmp=match(rownames(tt.down),genome$ensembl_gene_id)
	gene.down=data.frame(genome[tmp,c(1:2,4:7,13,16)],tt.down)
	write.cn(gene.down,paste0(name,".genelist.down.tsv"))

}

get_genename_4topgo=function(goID,GOdata,sig.gene,genome)
    {
	  library(topGO)
          gt=genesInTerm(GOdata,goID)
          tmp=(gt[[1]] %in% sig.gene)
          g=gt[[1]][tmp]          
          #gene_name = as.vector(genome$Associated.Gene.Name)[match(g,genome$Ensembl.Gene.ID)]
		  gene_name = as.vector(genome$external_gene_name)[match(g,genome$ensembl_gene_id)]
          return(paste(gene_name,collapse=","))
      } 
go.topgo=function(tt, name, score.cutoff, log2fc=0, mapping, genome)
{
	library(topGO)

	#tt=ttls[[ttname]]
	#gene.score=tt$PValue
	gene.score=tt$FDR
	names(gene.score)=rownames(tt)
		
	topdiff.up=function(allScore)
	{
		gene.score < score.cutoff & tt$logFC > log2fc
	}
	topdiff.down=function(allScore)
	{
		gene.score < score.cutoff & tt$logFC < -log2fc
	}
	
    test=function(GOdata)
        {
            des=description(GOdata)
            message(des)
            resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
			s = score(resultFisher)			
            allRes <- GenTable(GOdata, classicFisher = resultFisher,
                               orderBy = "classicFisher",
                               ranksOf = "classicFisher",topNodes = length(which(s<0.01)))
            allRes=subset(allRes,Significant > 0)
            sig.gene=sigGenes(GOdata) 
	
            gls=lapply(allRes$GO.ID,get_genename_4topgo,GOdata,sig.gene,genome)
            allRes$genes=unlist(gls)        
            write.cn(allRes,paste0("go.",des,".tsv"))
            save(GOdata,resultFisher,allRes,file=paste0("go.",des,".rda"))  
        }

	godata.up = new ("topGOdata", description = paste0(name,".up"), ontology = "BP", 
					allGenes = gene.score, geneSel = topdiff.up, nodeSize=10, 
					annot = annFUN.org, mapping = mapping, ID="Ensembl")
	test(godata.up)
	godata.down = new ("topGOdata", description = paste0(name,".down"), ontology = "BP", 
					allGenes = gene.score, geneSel = topdiff.down, nodeSize=10, 
					annot = annFUN.org, mapping = mapping, ID="Ensembl")    
    test(godata.down)
}

plotma_edger=function(tt,name,fdr=0.01, logfc=0)
    {
		pngname=paste0(name,".MA.png")
		png(pngname)
			plot(tt$logCPM,tt$logFC,pch=20,xlab="log2CPM",ylab="log2FC",main=name,col="grey")
			up=subset(tt,FDR< fdr & logFC> logfc)
			down=subset(tt,FDR< fdr & logFC< -logfc)
			points(up$logCPM,up$logFC,col="red",pch=20)
			points(down$logCPM,down$logFC,col="green",pch=20)
			legend("topright"
				,legend=c(paste0("FDR<",fdr),paste0("abs(logFC)>",round(logfc,2)),paste0("up=",nrow(up)),paste0("down=",nrow(down)))
				,col=c("white","white","red","green")
				,pch=20)
		dev.off()
    }

maplotSVG<-function(tt, output, fdr=0.05, log2fc=0,genome) 
{
library(htmlwidgets)
library(metricsgraphics)
library(RSVGTipsDevice)

id=rownames(tt)
ind=match(id,genome[,2])
tt$genename=genome[ind,1]
#ind=match(id,genome[,1])
#tt$genename=genome[ind,2]
tt$biotype=genename=genome[ind,13]

ind.up=which(tt$FDR<=fdr & tt$logFC>log2fc)
ind.down=which(tt$FDR<=fdr & tt$logFC< -log2fc)

imgname=paste0(output,".MA.svg")
f=function(i,cols)
{
	
    setSVGShapeToolTip(title=paste("Gene=",tt$genename[i], ": Biotype=",tt$biotype[i],": FDR=",round(tt$FDR[i],4), sep=""))
    points(tt$logCPM[i], tt$logFC[i], col=cols)
}

devSVGTips(imgname, toolTipMode=1, title=output, width =8, height = 8)  # this title is for the web page
    plot(tt$logCPM, tt$logFC, col="grey", xlim=range(tt$logCPM), ylim=range(tt$logFC), xlab="log2CPM", ylab="log2FC", main=output ) # empty plot frame
    abline(h=0)
	legend("topright"
		,legend=c(paste0("FDR<",fdr),paste0("abs(logFC)>",round(log2fc,2)),paste0("up=",length(ind.up)),paste0("down=",length(ind.down)))
		,col=c("white","white","red","green")
		,pch=1
		) 
				
	invisible(
		{
			sapply(ind.up,f,"red")
			sapply(ind.down,f,"green")
		}
	)
dev.off()


}


edger.chip=function(input,column,group,annpgr,output,norep=F,norm=F,fdr=0.01, logfc=0)
{
	library(edgeR)
	
	load(paste0(input,".chip.countls"))
	counts=count.table[,column]	
	group=as.factor(group)
	
	dge <- DGEList(counts, group=group)
	if (norm==T) dge <- calcNormFactors(dge) else
	{
		dge$samples$lib.size=size[column]
		dge <- calcNormFactors(dge,"none")
	}
	design <- model.matrix( ~ group, data=dge$samples )
	#design <- model.matrix( ~ group + batch, data=dge$samples )
	if (norep==F)
	{
	dge <- estimateGLMCommonDisp(dge, design)
	dge <- estimateGLMTrendedDisp(dge, design)
	dge <- estimateGLMTagwiseDisp(dge, design)
	fit <- glmFit(dge, design, prior.count=0)
	} else
	fit <- glmFit(dge, design, prior.count=0,dispersion=0.01 )
	lrTest <- glmLRT(fit, coef=2)
    #lrTest=lapply(list(2,3,4,5),glmLRT,glmfit=fit)
    #lrTest <- glmLRT(fit, contrast=c(1,-1,0,0,0))
    tt=lrTest$table
	tt$FDR=p.adjust(tt$PValue,method="BH")
										#tt=topTags(lrTest,Inf)
                                        #tt=lapply(lrTest,topTags,Inf)
                                        #x=lapply(tt,"[[",1)
	
		if (norep ==F)
		{
		png(paste0(output,".MDS.png"))
		plotMDS.DGEList(dge,col=as.numeric(group))
		dev.off()	
		}
		
#		fdr=0.01
#		logfc=0
#		png(paste0(output,".chip.MA.png"))
#		plot(tt$logCPM,tt$logFC,pch=20,xlab="log2CPM",ylab="log2FC",main=output,col="grey")
#			up=subset(tt,FDR< fdr & logFC> logfc)
#			down=subset(tt,FDR< fdr & logFC< -logfc)
#			points(up$logCPM,up$logFC,col="red",pch=20)
#			points(down$logCPM,down$logFC,col="green",pch=20)
#			legend("topright"
#				,legend=c(paste0("FDR<",fdr),paste0("abs(logFC)>",round(logfc,2)),paste0("up=",nrow(up)),paste0("down=",nrow(down)))
#				,col=c("white","white","red","green")
#				,pch=20)
#		dev.off()
	plotma_edger(tt,paste0(output,".chip"),fdr=fdr, logfc=logfc)	
			
	# out.edger=data.frame(annpgr[,c("seqnames","peak_start","peak_end","peak_width","gene.id","strand","genename","biotype","RefSeq.mRNA","peak_to_feature")],counts,count.norm[,column],tt)
	out.edger=data.frame(annpgr,counts,count.norm[,column],tt)
	save(dge,lrTest,tt,out.edger,file=paste0(output,".chip.edger"))
	write.cn(out.edger,paste0(output,".chip.edger.tab"))
	
}

if (FALSE)		   
{ # external example
group <- factor(c(1,1,2,2))
dge <- DGEList(counts, group=group)
dge <- calcNormFactors(dge)

## filter uniformative genes
m <- 1e6 * t(t(dge$counts) / dge$samples$lib.size)
ridx <- rowSums(m > 1) >= 2
dge <- dge[ridx,]

## comparison between groups
design <- model.matrix( ~ group)
#design <- model.matrix( ~ 0+group , data=dge$samples )
# colnames(design)=levels(dge$samples$group)
dge <- estimateCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

plotMDS.DGEList(dge)

## exact test
et <- exactTest(dge, pair=c("1","2"))
	topTags(et)
	detags <- rownames(topTags(et)$table)
	cpm(dge)[detags, order(dge$samples$group)]
	summary(de <- decideTestsDGE(et, p=1e-5))  # p is ?
	detags <- rownames(dge)[as.logical(de)]
	plotSmear(et, de.tags=detags)

## glm fit
fit <- glmFit(dge, design, dispersion=dge$common.dispersion)
lrTest <- glmLRT(fit, coef=2)
tt <- topTags(lrTest, Inf)
save(tt, file=file.path(dataDir, "tt.rda"))
}

if (FALSE)
{    
deseq.test=function(counts,group,output)
{
			library(DESeq)
			condition=as.factor(group)
           cds <- newCountDataSet(counts, condition)
           cds = estimateSizeFactors( cds )
		   #sizeFactors( cds )=size/mean(size)
		   cds = estimateDispersions( cds )
		   #plotDispEsts( cds )
		   res = nbinomTest( cds, levels(condition)[1],levels(condition)[2] )
		   plotMA(res)
#distance heatmap			
			cdsBlind= estimateDispersions( cds, method="blind" )
			vsd=varianceStabilizingTransformation( cdsBlind )
			dists = dist( t( exprs(vsd) ) )
			mat = as.matrix( dists )
			#rownames(mat) = colnames(mat) = with(pData(cdsBlind), paste(condition, batch, sep=" : "))
			rownames(mat) = colnames(mat) = with(pData(cdsBlind), paste(condition))
			library("RColorBrewer")
			hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
			library("gplots")
			heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))
			
			plotPCA(vsd)
          
#with batch		   
	design=data.frame(condition=group,batch=batch)
	cds=newCountDataSet(counts,design)
	cds = estimateSizeFactors( cds )
	#sizeFactors( cds )=size/mean(size)
	cds = estimateDispersions( cds )
	#plotDispEsts( cds )
	fit1=fitNbinomGLMs(cds,count ~ batch + condition)
	fit0=fitNbinomGLMs(cds,count ~ batch )
	pvalsGLM = nbinomGLMTest( fit1, fit0 )
	padjGLM = p.adjust( pvalsGLM, method="BH" )
	
	cdsBlind= estimateDispersions( cds, method="blind" )
	vsd=varianceStabilizingTransformation( cdsBlind )
	plotPCA(vsd, intgroup=c("condition", "batch"))
	
}


#pvalue of chip peak change, poisson p
chipdiff=function(countls,group)
{
#scale counts by size factor
	message("normalizing")
	counts=do.call(cbind,lapply(countls,"[[",1))
	size=unlist(lapply(countls,"[[",2))
	counts_norm=t(t(counts)/size)*1e6
	colnames(counts_norm)=paste("norm",colnames(counts_norm),sep=".")
	output=cbind(as.data.frame(tss.gr),counts,counts_norm)
	
	group=as.factor(group)
	g=as.list(levels(group))
	f.avg=function(x)
	{
	cols=which(group==x)
	if (length(cols)>1) average=rowMeans(counts_norm[,cols]) else
		average=counts_norm[,cols]
	return(average)
	}

	avg=lapply(g,f.avg)
	names(avg)=levels(group)
	tmp=do.call(cbind,avg)
	tmp1=log2(tmp[,-1]/tmp[,1])
	colnames(tmp1)=paste("log2fc",colnames(tmp1),sep=".")
	
	output=cbind(output,tmp,tmp1)
	
#pvalue
	counts_norm_size=t(t(counts)/size)*mean(size)
	wt=which(group==g[[1]])
	
	f.p=function(x)
	{
	cols=which(group==x)
	p1=ppois(rowMeans(counts_norm_size[,cols]),rowMeans(counts_norm_size[,wt]))
    p2=ppois(rowMeans(counts_norm_size[,cols]),rowMeans(counts_norm_size[,wt]), lower.tail = F)
	tmp=paste("log2fc",x,sep=".")
	col.fc=which(colnames(output)== tmp)
	p=ifelse(output[,col.fc]<0,p1,p2)	
		return(p)
	}

	p=lapply(as.list(levels(group)[-1]),f.p)
	names(p)=paste("pvalue",levels(group)[-1],sep=".")
	tmp=do.call(cbind,p)
	
	
#FDR
	#fdr=p*nrow(output)/rank(p)
	#fdr1=p1*nrow(output)/rank(p1)
	#output$FDR=ifelse(output$log2fc<0,fdr,fdr1)
	
	output=cbind(output,tmp)
	output
}

}
