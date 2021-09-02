args=commandArgs(trailingOnly=TRUE)

#Alias
rn<-rownames;
cn<-colnames;

#check args
if(length(args)!=5){
	stop("5 arguments needed: Rscript filter.R <expression_file> <output_folder> <number_of_replicates> <min genes> <min cov>\n",call.=FALSE)
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

cv<-function(x){ #Coefficient of variation
	return(sd(x)/mean(x));
}

##### specifics functions

#statistics on each sample
examineRNAseqSamples<-function(x, uncenter= FALSE){
	if(uncenter){
		x<-x-min(x)
		zero<-0
	}else{
		zero<-min(x)
	}
	mean<-apply(x,2,mean)
	sd<-apply(x,2,sd)
	count<-colSums(x)
	CV<-apply(x,2,cv)
	noGenEx<-rep(0,ncol(x))
	for(i in 1:ncol(x)) noGenEx[i]<-length(which(x[,i]>zero))

	return(data.frame(mean=mean, sd=sd,CV=CV,TotalGenEx=noGenEx,TotalCount=count))
}

filterMostExprimedGenesBySample <-function(data, numberOfGenes=nrow(data),minFreqOfGene=1,maxFreqOfGene=ncol(data),threshold=min(data)){
	#Fonction de filtrage des counts, on prend les x gènes les plus exprimés pour chache échantillon puis on effectue une jointure complète
	#puis on filtre selon la fréquence du gène
	#data: dataframe des counts
	#numberOfGenes: cutoff sur les rank d'expression des gènes, plus il est élevé moins le filtre est stringeant
	#minFreqOfGene : nombre de fois minimum où le gène revient dans la liste des x gènes les plus exprimés (où le gène > threshold)
    	#maxFreqOfGene : nombre de fois maximale où le gène est exprimé
	#threshold: l'expression du gène doit être plus grande que ce paramètre pour que la fréquence soit comptée comme 1
	mostExprimedGenes<-list()
	for(i in 1:ncol(data)){
		col<- data[,i]
		names(col)<-rownames(data)
		col<-col[which(col>threshold)];
		mostExprimedGenes[[i]]<-names(sort(col,decreasing = T))[1:min(numberOfGenes,length(col))]
	}
	rm(col)

	freqTable<-summary(as.factor(unlist(mostExprimedGenes)),maxsum=nrow(data))
	mostExprimedgenesVector<-names(freqTable[which(freqTable>=minFreqOfGene & freqTable<=maxFreqOfGene)])
	return(data[mostExprimedgenesVector,])
}


# start
expressionData=args[1]
outputDir=args[2]
nbRep=as.integer(args[3])
minimumGeneExpressed=as.integer(args[4])
minimumTotalCount=as.integer(args[5])

exprDat<-lire(expressionData)

absSamples<-examineRNAseqSamples(exprDat)
ecrire(absSamples,paste(outputDir,"/SamplesAbstract.tsv",sep=""))

sample2keep<-rn(absSamples)[absSamples$TotalGenEx>minimumGeneExpressed & absSamples$TotalCount>minimumTotalCount]
write.table(sample2keep,paste(outputDir,"/samplesToKeep.txt",sep=""),col.names=F,row.names=F,quote=F)
sample2rm<-cn(exprDat)[!cn(exprDat)%in%sample2keep]; 

if(length(sample2rm)>0){
	warning("these samples don't pass the quality control, they will be removed from analysis: ");warning(paste(sample2rm," ",sep=" "))
}

exprDat<-exprDat[sample2keep]

exprDat<-filterMostExprimedGenesBySample(exprDat,minFreqOfGene=nbRep) # Gene filtering

ecrire(exprDat,paste(outputDir,"/exprFiltered.tsv",sep=""))
ecrire(absSamples[sample2rm,c("TotalGenEx","TotalCount")],paste(outputDir,"/samplesToRemove.txt",sep=""))


