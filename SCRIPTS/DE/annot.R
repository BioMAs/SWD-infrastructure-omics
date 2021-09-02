args=commandArgs(trailingOnly=TRUE)
outputDir=args[1]
cond1=args[2]
cond2=args[3]
org=args[4]
corresFile=args[5]

library(GO.db)
library(GSEABase)
library(clusterProfiler)
library(DOSE)
library(fgsea)

makeGeneList<-function(data,org){
	pAdj=data[,"padj"]
	InvPAdj=1-pAdj
	data=data.frame(data, InvPAdj)
	logC=data[,"log2FoldChange"]
	logC[logC<0] <- -1
	logC[logC>=0] <- 1
	data$InvPAdj=data$InvPAdj*logC
	geneID=bitr(data$Name, fromType="SYMBOL",toType=c("ENTREZID"),OrgDb=org)
	doublons=which(duplicated(geneID$SYMBOL))
	if(length(doublons)>0){geneID=geneID[-doublons,]}
	data=data[is.element(data$Name, geneID$SYMBOL),]
	geneList=data[,dim(data)[2]]
	names(geneList)=as.character(geneID$ENTREZID)
	geneList=sort(geneList, decreasing=TRUE)		
	return(geneList)
}

corresIDorg=read.table(corresFile,header=T,row.names=1)
orgGO=as.character(corresIDorg[org,1])
orgKegg=as.character(corresIDorg[org,2])

comp=paste(cond1,cond2,sep="__vs__")

data=read.csv(paste(paste(outputDir,comp,sep="/"),"DEseqResFiltered.tsv",sep="/"), header=TRUE, sep="\t") 

#Sys.setenv(http_proxy="http://cache.ha.univ-nantes.fr:3128")
#Sys.setenv(https_proxy="https://cache.ha.univ-nantes.fr:3128")

if(dim(data)[1]>10){
	geneList=makeGeneList(data,orgGO)
	gene=names(geneList)
	ego=enrichGO(gene=gene,OrgDb=orgGO,ont="ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05,readable      = TRUE)
	png(paste(paste(outputDir,comp,sep="/"),"dotplotGO.png",sep="/"),width=1000,height=600)
	print(dotplot(ego))
	dev.off()
	write.table(ego,paste(paste(outputDir,comp,sep="/"),"annotGo.txt",sep="/"),sep="\t",row.names=F,quote=F)

	ekegg=enrichKEGG(gene, organism = as.character(orgKegg), keyType = "kegg", pvalueCutoff = 0.05,pAdjustMethod = "BH", minGSSize = 10, maxGSSize = 500,qvalueCutoff = 0.2, use_internal_data = FALSE)
	eKeggSymbol=setReadable(ekegg,orgGO,keyType="ENTREZID")
	png(paste(paste(outputDir,comp,sep="/"),"dotplotKEGG.png",sep="/"),width=1000,height=600)
	print(dotplot(eKeggSymbol))
	dev.off()
	write.table(eKeggSymbol,paste(paste(outputDir,comp,sep="/"),"annotKegg.txt",sep="/"),sep="\t",row.names=F,quote=F)

}


dataTot=read.csv(paste(paste(outputDir,comp,sep="/"),"DEseqRes.tsv",sep="/"), header=TRUE, sep="\t")
geneListTot=makeGeneList(dataTot,orgGO)

ego2Symbol <- tryCatch(
  {
    ego2 <- gseGO(geneList=geneListTot,OrgDb=orgGO,ont="ALL",nPerm=1000,minGSSize=20,maxGSSize=500,pvalueCutoff=1,verbose=FALSE)
    setReadable(ego2,orgGO,keyType="ENTREZID")
  },
  error = function(e){
    matrix(nrow = 1, ncol = 12, dimnames=list('',c('ONTOLOGY','ID','Description','setSize','enrichmentScore','NES','pvalue','p.adjust','qvalues','rank','leading_edge','core_enrichment')) )
  }
)


kk2Symbol <- tryCatch(
  {
    kk2 <- gseKEGG(geneList=geneListTot,organism=as.character(orgKegg),nPerm=1000,minGSSize=2,pvalueCutoff=1,verbose=FALSE)
    setReadable(kk2,orgGO,keyType="ENTREZID")
  },
  error = function(e){
    matrix(nrow = 1, ncol = 11, dimnames=list('',c('ID','Description','setSize','enrichmentScore','NES','pvalue','p.adjust','qvalues','rank','leading_edge','core_enrichment')) )
  }
)


write.table(ego2Symbol,paste(paste(outputDir,comp,sep="/"),"gseGo.txt",sep="/"),sep="\t",row.names=F,quote=F)
write.table(kk2Symbol,paste(paste(outputDir,comp,sep="/"),"gseKegg.txt",sep="/"),sep="\t",row.names=F,quote=F)

