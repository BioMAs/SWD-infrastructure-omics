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

##### specifics functions

UMI2UPM<-function(data){ #Normalisation UPM
        data.UPM <- sweep(data, 2, colSums(data),`/`)
        data.UPM <-data.UPM * 1000000
        return(data.UPM)
}


expressionData=args[1]
outputDir=args[2]
exprDat<-lire(expressionData)
exprDat.UPM=UMI2UPM(exprDat)
ecrire(exprDat.UPM,paste(outputDir,"/exprDatUPM.tsv",sep=""))

