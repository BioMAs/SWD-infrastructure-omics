library(DESeq2)
library(limma) #removeBatchEffect

args=commandArgs(trailingOnly=TRUE)

if(length(args)==3){
	batchCol=NULL
}else{
	batchCol=args[4]
}

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

#start
expressionData=args[1]
sampleTable=args[2]
outputDir=args[3]
condCol="Condition"

exprDat<-lire(expressionData)
sampleAnnot<-lire(sampleTable)

batch<-TRUE; if(is.null(batchCol)) batch<-FALSE

sampleAnnot[,condCol]<-as.factor(as.character(sampleAnnot[,condCol]))
exprDat<-exprDat[,c(rownames(sampleAnnot))]

if(length(unique(sampleAnnot[,condCol]))>1){
	formulaChar<-paste0("~",condCol)
}else{
	formulaChar<-"~1"
}


if(batch){
	sampleAnnot[,batchCol]<-as.factor(as.character(sampleAnnot[,batchCol]))
	formulaChar<-paste0(formulaChar,"+",batchCol)
}

dds <- DESeqDataSetFromMatrix(countData = exprDat,colData = sampleAnnot,design = formula(formulaChar)) 
dds <- DESeq(dds)

save(dds,file=paste(outputDir,"/dds.RData",sep=""))

exprDatN<-counts(dds,normalize=TRUE)

if(ncol(exprDat)<=20){
	exprDatT<-assay(rlog(dds))
}else{
	exprDatT<-assay(vst(dds))
}

ecrire(exprDatN,paste(outputDir,"/exprNormalized.tsv",sep=""))
ecrire(exprDatT,paste(outputDir,"/exprTransformed.tsv",sep=""))


if(batch){
	exprDatNa<-removeBatchEffect(exprDatN,batch = sampleAnnot[,batchCol],design = model.matrix(data=sampleAnnot,formula(paste0("~",condCol))))
	exprDatTa<-removeBatchEffect(exprDatT,batch = sampleAnnot[,batchCol],design = model.matrix(data=sampleAnnot,formula(paste0("~",condCol))))
	ecrire(exprDatNa,paste(outputDir,"/exprNormalizedAdjusted.tsv",sep=""))
	ecrire(exprDatTa,paste(outputDir,"/exprTransformedAdjusted.tsv",sep=""))
}

