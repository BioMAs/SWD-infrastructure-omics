args=commandArgs(T)
input=args[1]
output=args[2]
data=read.table(input,sep="\t", header=T, row.names=1, check.names=FALSE)

assigned=as.numeric(data[1,])
aligned=as.numeric(data[2,])
refseqTot=as.numeric(data[3,])
refSeqUMI=as.numeric(data[4,])

#mat=cbind(refSeqUMI,refseqTot,aligned,assigned)



pdf(output,width=10,height=6)
par(mar=c(10.1,4.1,4.1,6.1), xpd=TRUE)
barplot(assigned,col=c("black"),names.arg=colnames(data),cex.names=0.5,las=2,main="reads repartition",ylim=c(0,max(assigned)+500000),axes=F)
axis(2,seq(0,max(assigned)+500000,500000),format(seq(0,max(assigned)+500000,500000),big.mark=" ", scientific=F),las=2,cex.axis=0.5)
barplot(aligned,col=c("red"),add=T,yaxt="n")
barplot(refseqTot,col=c("green"),add=T,yaxt="n")
barplot(refSeqUMI,col=c("blue"),add=T,yaxt="n")
legend("right",inset=c(-0.1,0),legend=c("assigned","aligned","refseq Tot","refseq UMI"),col=c("black","red","green","blue"),lty=1,lwd=2,cex=0.7,bty="n")
dev.off()




