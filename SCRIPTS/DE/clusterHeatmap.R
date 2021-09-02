library(ggplot2) 
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
        return(gpar(fontsize=1/n*400))
}

corrDist<-function(x) return(as.dist((1 - cor(Matrix::t(x)))/2))

#start
exprDatTFile=args[1]
sampleTable=args[2]
sampleAbstract=args[3]
outputDir=args[4]
condCol="Condition"
absSamples=lire(sampleAbstract)

bootstrap<-FALSE #bootstrap unsupervised clustering ? bootstrap increase a lot computing time and may generate errors if N is too small
nboot=30  #if bootsrap, number of bootstrap replications

sampleAnnot<-lire(sampleTable)

if(length(args)==4){
        batchCol=NULL
	exprDatT=lire(exprDatTFile)
	annot<-sampleAnnot[condCol]	
}else{
        batchCol=args[3]
	exprDatT=lire(paste(outputDir,"/exprTransformedAdjusted.tsv",sep=""))
	annot<-sampleAnnot[,c(condCol,batchCol)]
}

corSample<-cor(exprDatT,method = "pearson")
clustSamples<-unsupervisedClustering(corSample,bootstrap = bootstrap,nboot=nboot,method.dist = "euclidean")

colTopAnnot<-vector("list", ncol(annot))
names(colTopAnnot)<-colnames(annot)
colFun<-c(ggplotColours,rainbow)
i<-1
for(col in colnames(annot)){
	colTopAnnot[[col]]<-colFun[[i]](nlevels(annot[,col]))
	names(colTopAnnot[[col]])<-levels(annot[,col])
	i<-i+1
	if(i==4) i<-1
}

ha<-HeatmapAnnotation(df = annot, col = colTopAnnot)

Ht<-Heatmap(matrix = corSample, cluster_rows = clustSamples, cluster_columns = clustSamples, top_annotation = ha,name="Pearson\ncorrelation",row_names_gp = autoGparFontSizeMatrix(ncol(corSample)),column_names_gp = autoGparFontSizeMatrix(ncol(corSample)),row_dend_reorder = FALSE, column_dend_reorder = FALSE,col=heat.colors(100))

Htc<-Heatmap(absSamples[clustSamples$labels,"TotalCount"], name = "Total\nexpression",col=colorRamp2(c(0,max(absSamples$TotalCount)),c("black","green")), show_row_names = FALSE, width = unit(5, "mm"))

pdf(file = paste(outputDir,"/HeatmapCorPearson.pdf",sep=""),width=10,height=9)

dev.set(1)
png(paste(outputDir,"/HeatmapCorPearson.png",sep=""),width=25,height=22,units="cm",res=600,pointsize=3)
print(Ht+Htc)
dev.off()

print(Ht+Htc)

if(!is.null(batchCol)){
	#il y a eu une correction batch
	exprDat=lire(paste(outputDir,"/exprTransformed.tsv",sep=""))
	corSampleUnadjust<-cor(exprDat,method = "pearson")
	clustSamplesUnadjust<-unsupervisedClustering(corSampleUnadjust,nboot=nboot,bootstrap = bootstrap,method.dist = "euclidean")

	Htna<-Heatmap(matrix = corSampleUnadjust, cluster_rows = clustSamplesUnadjust, cluster_columns = clustSamplesUnadjust, top_annotation = ha,name="Pearson\ncorrelation\nnon adjusted",row_names_gp = autoGparFontSizeMatrix(ncol(corSampleUnadjust)),column_names_gp = autoGparFontSizeMatrix(ncol(corSampleUnadjust)),row_dend_reorder = FALSE, column_dend_reorder = FALSE,col=heat.colors(100))
	print(Htna+Htc)

	dev.set(1)
	png(paste(outputDir,"/HeatmapCorPearson.png",sep=""),width=25,height=22,units="cm",res=600,pointsize=3)
	print(Htna+Htc)
	dev.off()
}

dev.off()

save(ha,batchCol,file=paste(outputDir,"/heatmap.RData",sep=""))
