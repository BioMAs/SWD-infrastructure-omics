args=commandArgs(trailingOnly=TRUE)
org=args[1]
corresFile=args[2]

corresIDorg=read.table(corresFile,header=T,row.names=1)
orgGO=as.character(corresIDorg[org,1])

if((!require(orgGO))){
	requireNamespace("BiocManager", quietly = TRUE)
	BiocManager::install(orgGO, updates = TRUE)
}
