

pheno <- data.frame(ID=paste("P",1:200,sep=""), STATUS=c(rep("case",100),rep("control",100)))
require(combinat)

netDir <- "~/Documents/Software/netDx/netDx/inst/extdata/example_nets"
dir.create(netDir)

# big net all cases
x <- t(combn(pheno$ID[1:10],2))
x <- data.frame(x[,1],x[,2],1)
write.table(x,file=sprintf("%s/BIG_CASE.txt",netDir),sep="\t",col=F,row=F,
			quote=F)

# big net all controls
x <- t(combn(pheno$ID[101:110],2))
x <- data.frame(x[,1],x[,2],1)
write.table(x,file=sprintf("%s/BIG_CONTROL.txt",netDir),sep="\t",col=F,row=F,
			quote=F)

# half cases and half controls
x <- t(combn(pheno$ID[c(1:5,101:105)],2))
x <- data.frame(x[,1],x[,2],1)
write.table(x,file=sprintf("%s/BOTH_EQUAL.txt",netDir),sep="\t",col=F,row=F,
			quote=F)

# 3/4 cases
x <- t(combn(pheno$ID[c(1:7,101:103)],2))
x <- data.frame(x[,1],x[,2],1)
write.table(x,file=sprintf("%s/MOSTLY_CASE.txt",netDir),sep="\t",col=F,row=F,
			quote=F)

# small net all cases
x <- t(combn(pheno$ID[1:3],2))
x <- data.frame(x[,1],x[,2],1)
write.table(x,file=sprintf("%s/SMALL_CASE.txt",netDir),sep="\t",col=F,row=F,
			quote=F)

# small net all controls
x <- t(combn(pheno$ID[101:103],2))
x <- data.frame(x[,1],x[,2],1)
write.table(x,file=sprintf("%s/SMALL_CONTROL.txt",netDir),sep="\t",col=F,row=F,
			quote=F)
