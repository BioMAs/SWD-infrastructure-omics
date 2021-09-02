library(ggplot2) #pour fonction acp2d et plotDistrib
library(grid) #pour fonction plotDistrib

args=commandArgs(trailingOnly=TRUE)

exprDatTFile=args[1]
sampleTable=args[2]
outputDir=args[3]

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

# specific function
ACP<-function(d,transpose=T,scale=F,center=T) {
        if(transpose) d<-t(d);
        resacp <-prcomp(x = d,retx = T,center = center,scale = scale);
        resacp$n.obs<-dim(d)[1];
        resacp$percentVar<- resacp$sdev^2 / sum( resacp$sdev^2 )
        return(resacp);
}

acp2d<-function(pca, comp=1:2,group=NULL, plotVars = FALSE, pointSize=2, plotText=FALSE,fixedCoord=FALSE,main=NULL,ellipse=FALSE,color=NULL){
        if(!require("ggplot2")) stop("You must install ggplot2");
        if(length(comp)!=2) stop("You must give a vector of 2 integer for comp parameter");
        percentVar <- pca$percentVar
        functPlot<-ifelse(plotText,geom_text,geom_point)
        coor=ifelse(plotVars,"rotation","x")

        if(is.null(group)){
                d <- data.frame(PC1=pca[[coor]][,comp[1]], PC2=pca[[coor]][,comp[2]]);
                graph<-ggplot(data=d, mapping = aes(x=PC1, y=PC2, label = rownames(d)))
        }else{
                d <- data.frame(PC1=pca[[coor]][,comp[1]], PC2=pca[[coor]][,comp[2]], group=group);
                graph<-ggplot(data=d, mapping = aes(x=PC1, y=PC2,colour=group, label = rownames(d)))
        }
        graph<-graph+functPlot(size=pointSize)+
                xlab(paste0("PC",comp[1],": ",round(percentVar[comp[1]] * 100),"% variance")) +
                ylab(paste0("PC",comp[2],": ",round(percentVar[comp[2]] * 100),"% variance"))
        if(fixedCoord)graph <- graph + coord_fixed(ratio=percentVar[comp[2]]/percentVar[comp[1]])
        if(ellipse) graph<-graph+stat_ellipse()
        if(!is.null(main)) graph <- graph + ggtitle(main)
        if(!is.null(color)) graph <- graph + scale_color_manual(values=color)

        print(graph)
}

plotDistrib<-function(data,type="boxplot",conditions=NULL,main=NULL,conditionName="Batch"){
  require(grid)
  require(ggplot2)
  #data : matrix of expression data
  #type: 'violin' or 'boxplot' ?
  #conditions : vector of factors to color plot
  if(!type%in%c("violin","boxplot")) stop("type must be 'violin', 'hist' or 'boxplot'")

  vectorDat<-as.vector(as.matrix(data))
  tabGraph<-data.frame(val=vectorDat,sample=rep(colnames(data),each=nrow(data)))
  if(!is.null(conditions)){
    tabGraph[,conditionName]=rep(conditions,each=nrow(data))
  }

  tabGraph$sample<-factor(tabGraph$sample,levels=colnames(data))


  if(is.null(conditions)){
    graph<-ggplot(data = tabGraph,mapping = aes(sample,val))
  }else{
    graph<-ggplot(data = tabGraph,mapping = aes_string("sample","val",color=conditionName))
  }

    graph<-graph+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=.3))
  if(type=="violin") graph<-graph+geom_violin()
  if(type=="boxplot") graph<-graph+geom_boxplot()
  if(! is.null(main)) graph<-graph+ggtitle(main)

  print(graph)
}

plotAllACP<-function(compo=4,acp,sampleAnnot,condCol){
	for(i in 1:(compo-1)){
  		for(j in (i+1):compo){
    			if(i == 1 & j ==2){
				acp2d(acp,group = sampleAnnot[,condCol],plotText = TRUE,pointSize = 4,comp = c(i,j),main="PCA",fixedCoord = F)
			}
    			acp2d(acp,group = sampleAnnot[,condCol],pointSize = 2,comp = c(i,j),main="PCA",fixedCoord = F)
    			#acp2d(acp,pointSize = 2,comp = c(i,j),plotVars = TRUE, plotText = TRUE,fixedCoord = F)
  		}
	}
}

#start
#exprDatT=lire(paste(outputDir,"/exprTransformed.tsv",sep=""))
exprDatT=lire(exprDatTFile)
condCol="Condition"
sampleAnnot<-lire(sampleTable)
acpT<-ACP(exprDatT)
compo=length(which(acpT$percentVar>0.05))+1
exprDatT<-exprDatT[,c(rownames(sampleAnnot))]

pdf(file = paste(outputDir,"/NormAndPCA.pdf",sep=""),width=10,height=10)

if(!is.null(batchCol)){
        #il y a eu une correction batch
	exprDatTa=lire(paste(outputDir,"/exprTransformedAdjusted.tsv",sep=""))
	acpTa<-ACP(exprDatTa)
	plotDistrib(exprDatT,main="Boxplots of normalized/transformed data",conditions = sampleAnnot[,batchCol])
	plotDistrib(exprDatTa,main="Boxplots of normalized/transformed/adjusted data",conditions = sampleAnnot[,batchCol])
	barplot(acpTa$percentVar*100,names.arg = round(acpTa$percentVar*100,2),main = "Contribution of each componant in PCA")
	
	acp2d(acpT,comp=c(1,2),plotText = TRUE,group = sampleAnnot[,batchCol],main="PCA on normalized/transformed data")
	#acp2d(acpTa,comp=c(1,2),plotText = TRUE,group = sampleAnnot[,condCol],main="PCA on normalized/transformed/adjusted data")
	acp2d(acpTa,comp=c(1,2),plotText = TRUE,group = sampleAnnot[,batchCol],main="PCA on normalized/transformed/adjusted data")

	dev.set(1)
	png(paste(outputDir,"/PCA.png",sep=""),width=15,height=10,units="cm",res=300,pointsize=6)
	acp2d(acpTa,comp=c(1,2),group = sampleAnnot[,batchCol],main="PCA on normalized/transformed/adjusted data")
	dev.off()

	plotAllACP(compo,acpTa,sampleAnnot,condCol)

}else{
	#pas de correction batch
	plotDistrib(exprDatT,main="Boxplots of normalized/transformed data")
	#acp2d(acpT,comp=c(1,2),plotText = TRUE,group = sampleAnnot[,condCol],main="PCA on normalized/transformed data")
	barplot(acpT$percentVar*100,names.arg = round(acpT$percentVar*100,2),main = "Contribution of each componant in PCA")
	plotAllACP(compo,acpT,sampleAnnot,condCol)

	dev.set(1)
        png(paste(outputDir,"/PCA.png",sep=""),width=15,height=10,units="cm",res=300,pointsize=6)
        acp2d(acpT,comp=c(1,2),group = sampleAnnot[,condCol],main="PCA on normalized/transformed data")
        dev.off()

}

dev.off()


