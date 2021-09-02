library(fdrtool)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(pvclust)
library(circlize)


args=commandArgs(trailingOnly=TRUE)

# usual functions
lire<-function(x, character=FALSE){
        if(character){
                d<-read.table(file = x,sep = "\t",header=T,row.names = 1,colClasses = "character",quote="",check.names=FALSE)
        }else{
                d<-read.table(file = x,sep = "\t",header=T,row.names = 1,quote="",check.names=FALSE)
        }
        return(d)
}

ecrire<-function(x,file="default.tsv",headRow="Name"){
        options(warn=-1) #Supress unecessary warning about append
        write.table(x = paste0(headRow,"\t"),file = file,sep = "\t",eol="",quote=F,row.names=F,col.names=F)
        write.table(x=x,file=file,sep="\t", row.names = T, col.names = T, quote = FALSE,append=T)
        options(warn=0)
}

# specific function
unsupervisedClustering<-function(x,transpose=TRUE,method.dist="pearson",method.hclust="ward.D2",bootstrap=FALSE,nboot=10){
        if(transpose) x<-t(x)
        if(bootstrap){
                require(pvclust)
                resClust<-pvclust(t(x),nboot=nboot,method.hclust = method.hclust,parallel = TRUE,method.dist = method.dist)$hclust
        }else{
                if(method.dist=="pearson"){
                        resDist<-corrDist(x)
                }else{
                        resDist<-dist(x, method = method.dist)
                }
                resClust<-stats::hclust(resDist,method = method.hclust)
        }
        return(resClust)
}

ggplotColours <- function(n = 6, h = c(0, 360) + 15){
    if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

autoGparFontSizeMatrix<-function(n){ #Calcule automatiquement la taille de police selon le nombre de colonnes ou lignes (empirique)
        n=max(n,50)
        n=min(n,1000)
        return(gpar(fontsize=1/n*600))
}

rowScale<-function(data,center=TRUE,scaled=FALSE){
        data<-t(data)
        data<-t(scale(data,center=center,scale=scaled))
        return(data)
}

corrDist<-function(x) return(as.dist((1 - cor(Matrix::t(x)))/2))

#start
outputDir=args[1]
cond1=args[2]
cond2=args[3]
condCol="Condition"
load(paste(outputDir,"/dds.RData",sep=""))
load(paste(outputDir,"/heatmap.RData",sep=""))
logFCthreshold=as.numeric(args[4])
AdjPValthreshold=0.05
GenesInFig=50
bootstrap=FALSE
nboot=30

sampleAnnot<-lire(paste(outputDir,"/sampletable.tsv",sep=""))
comp<-paste(cond1,cond2,sep="__vs__")

samples<-rownames(sampleAnnot)

res=results(dds,contrast=c(condCol,cond1,cond2),independentFiltering=T)

if(file.exists(paste(outputDir,"/exprNormalizedAdjusted.tsv",sep=""))){
	exprDat=lire(paste(outputDir,"/exprNormalizedAdjusted.tsv",sep=""))
	exprDatT=lire(paste(outputDir,"/exprTransformedAdjusted.tsv",sep=""))
	batch=TRUE
}else{
	exprDat=lire(paste(outputDir,"/exprNormalized.tsv",sep=""))
	batch=FALSE
	exprDatT=lire(paste(outputDir,"/exprTransformed.tsv",sep=""))
}
res$meanInComp<-rowMeans(exprDat[,sampleAnnot[,condCol]%in% c(cond1,cond2)])

DE=data.frame(res)

DE.sel=list()
DE.sel$up=DE[which(DE$padj<AdjPValthreshold & DE$log2FoldChange > logFCthreshold),]
DE.sel$down=DE[which(DE$padj<AdjPValthreshold & DE$log2FoldChange < -logFCthreshold),]
DE.sel$isDE=rbind(DE.sel$up,DE.sel$down)
DE.sel$notDE=DE[setdiff(rownames(DE),rownames(DE.sel$isDE)),]
DE$DE="NONE"
DE[rownames(DE.sel$up),"DE"]="UP"
DE[rownames(DE.sel$down),"DE"]="DOWN"
DE$DE=factor(DE$DE,levels=c("DOWN","NONE","UP"))

ecrire(DE,paste(paste(outputDir,comp,sep="/"),"DEseqRes.tsv",sep="/"))
ecrire(rbind(DE.sel$up,DE.sel$down),paste(paste(outputDir,comp,sep="/"),"DEseqResFiltered.tsv",sep="/"))

#Volcano-plot
data <- cbind(DE, color = DE$DE)
levels(data$color)[levels(data$color) == "DOWN"] <- "green"
levels(data$color)[levels(data$color) == "NONE"] <- "grey75"
levels(data$color)[levels(data$color) == "UP"] <- "red"
data <- cbind(data, label = ifelse(data$DE != "NONE", rownames(data), ""))

# removing NA for setting correct plot limits
data <- data[!is.na(data$padj),]
# limiting number of displayed genes
data <- data[order(abs(data$padj)),]
data[GenesInFig+1:nrow(data),]$label <- ""
dataLabels <- subset(data, label != "")

png(paste(paste(outputDir,comp,sep="/"),"Volcano-plot.png",sep="/"),width=10,height=10,units="cm",res=200)

p <- ggplot(data=data, aes(x=log2FoldChange, y=-log10(padj))) +
    geom_point(size = 0.5, color = data$color) +
    geom_hline(aes(yintercept=-log10(AdjPValthreshold))) +
    geom_vline(aes(xintercept=logFCthreshold)) +
    geom_vline(aes(xintercept=-logFCthreshold)) +
    ylab("-log10(DEseq padj)") +
    xlab("log2(Fold-Change)") +
    ggtitle(paste0("Volcano plot of comparison\n", cond1," (pos FC) VS ",cond2," (neg FC)\nBenjamini & Hochberg method (",nrow(DE.sel$isDE),"/",nrow(DE)," DE genes)")) +
    theme_bw() +
    theme(plot.title = element_text(size=10),
          axis.title = element_text(size=10))
# put labels when existing
if (nrow(dataLabels)>0) {
	p <- p + geom_text_repel(
	  data = dataLabels,
	  aes(label = label, fontface="bold.italic"),
	  size = 2,
	  point.padding = unit(0.1, "lines")
	)
}
print(p)
dev.off()

#MA-plot
data<-data.frame(DE[res$meanInComp>0.5,])
data$NAME<-rownames(data)
data<-data[order(abs(data$pval)),]
dataDE<-data[1:GenesInFig,]
data$DE=factor(data$DE,levels=c("DOWN","NONE","UP"))
myColors <- c("green","grey75","red")
names(myColors) <- levels(data$DE)

png(paste(paste(outputDir,comp,sep="/"),"MA-plot.png",sep="/"),width=10,height=10,units="cm",res=200)
ggplot(data=data,mapping = aes(x=meanInComp,y=log2FoldChange,colour=DE))+scale_colour_manual(name = "DE",values = myColors,guide=FALSE)+geom_point(size=0.2)+scale_x_log10()+
       geom_text_repel(data = dataDE,mapping = aes(x=meanInComp,fontface="bold.italic",
       y=log2FoldChange,label=NAME),inherit.aes = FALSE,max.iter = 8000,size=2)+
       theme(panel.background = element_rect(fill = NA,colour="black"),
       panel.grid.major = element_line(colour = NA),
       panel.grid.minor = element_line(colour = NA),plot.title = element_text(size=10),axis.title=element_text(size=10))+
       ggtitle(paste0("MA-plot for ",cond1," vs ",cond2))+
       xlab("Mean expression")+ylab("log2(Fold-Change)")
dev.off()

#DE genes clustering
if(nrow(DE.sel$isDE)>=10){

	DEgenes.names=rownames(DE.sel$isDE)
	colFun<-c(ggplotColours,rainbow)

	sampleHt<-colnames(exprDatT)
        haByComp<-ha
        exprDE=exprDatT[DEgenes.names,sampleHt,drop=FALSE]
        exprDE.scaled=rowScale(exprDE,center=T,scaled=T)

        bootTemp=bootstrap
        if(nrow(exprDE.scaled>10)){bootTemp=FALSE}

        hclustGeneDE<-unsupervisedClustering(exprDE.scaled,transpose = F,nboot=nboot,bootstrap = bootTemp)
        hclustSampleDE<-unsupervisedClustering(exprDE.scaled,transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")

        rowScaledExpr<-rowScale(exprDatT,center=TRUE,scale=TRUE)
        quantile.expr<-quantile(unlist(rowScaledExpr),seq(0,1,.01))
        colHA=colorRamp2(c(quantile.expr[2],0,quantile.expr[100]),c("blue","white","red"))

	pdf(paste(paste(outputDir,comp,sep="/"),"clustDEgene.pdf",sep="/"),width=10,height=10)
	p <- Heatmap(exprDE.scaled,top_annotation = haByComp, row_names_gp = autoGparFontSizeMatrix(nrow(exprDE.scaled)),
                cluster_rows = hclustGeneDE,col = colHA, cluster_columns = hclustSampleDE,name=comp,column_names_gp = autoGparFontSizeMatrix(ncol(exprDE.scaled)) )
    print(p)

	if(length(samples[sampleAnnot[,condCol]%in% c(cond1,cond2)])<length(sampleHt)){
		sampleHtSelect=samples[sampleAnnot[,condCol]%in% c(cond1,cond2)]
		if(batch){
			tempAnnot<-droplevels(sampleAnnot[sampleHtSelect,c(condCol,batchCol)]) 
		}else{
      			tempAnnot<-droplevels(sampleAnnot[sampleHtSelect,condCol,drop=FALSE])
		}
		colTopAnnotTemp<-vector("list", ncol(tempAnnot))
    		names(colTopAnnotTemp)<-colnames(tempAnnot)
    		i<-1
    		for(col in colnames(tempAnnot)){
      			colTopAnnotTemp[[col]]<-colFun[[i]](nlevels(tempAnnot[,col]))
      			names(colTopAnnotTemp[[col]])<-levels(tempAnnot[,col])
      			i<-i+1
      			if(i==4) i<-1
    		}
   		haByCompSelect<-HeatmapAnnotation(df = tempAnnot,col = colTopAnnotTemp)
		exprDESelect=exprDatT[DEgenes.names,sampleHtSelect,drop=FALSE]
		exprDESelect.scaled=rowScale(exprDESelect,center=T,scaled=T)
		exprDESelect.scaled<-na.omit(exprDESelect.scaled)

		hclustGeneDESelect<-unsupervisedClustering(exprDESelect.scaled,transpose = F,nboot=nboot,bootstrap = bootTemp)
        	hclustSampleDESelect<-unsupervisedClustering(exprDESelect.scaled,transpose = T,nboot=nboot,bootstrap = bootTemp,method.dist = "euclidean")
		
		png(paste(paste(outputDir,comp,sep="/"),"clustDEgene.png",sep="/"),width=25,height=22,units="cm",res=600,pointsize=3)
		print(Heatmap(exprDESelect.scaled,top_annotation = haByCompSelect, row_names_gp = autoGparFontSizeMatrix(nrow(exprDESelect.scaled)),
		cluster_rows = hclustGeneDESelect,col = colHA, cluster_columns = hclustSampleDESelect,name=comp,column_names_gp = autoGparFontSizeMatrix(ncol(exprDESelect.scaled)) ))	
		dev.off()

		print(Heatmap(exprDESelect.scaled,top_annotation = haByCompSelect, row_names_gp = autoGparFontSizeMatrix(nrow(exprDESelect.scaled)),
                cluster_rows = hclustGeneDESelect,col = colHA, cluster_columns = hclustSampleDESelect,name=comp,column_names_gp = autoGparFontSizeMatrix(ncol(exprDESelect.scaled)) ))

  	}
  	else {	
		png(paste(paste(outputDir,comp,sep="/"),"clustDEgene.png",sep="/"),width=25,height=22,units="cm",res=600,pointsize=3)
		print(p)
		dev.off()
  	}
	dev.off()

}else{
	png(paste(paste(outputDir,comp,sep="/"),"clustDEgene.png",sep="/"))
	plot.new()
	text(0.5,0.5,"Less than 10 differential expressed gene for this comparison,\n clustering isn't relevant")
	dev.off()
}
