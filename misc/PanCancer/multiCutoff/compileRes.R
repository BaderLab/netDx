#' compile best result from all tumours
rm(list=ls())

.getBest <- function(d) {
	d <- do.call("rbind",d)
	x <- apply(d,2,which.max)
	val <- apply(d,2,max)
	nm <- rownames(d)[x]
	z <- cbind(val,nm,names(x))
	z
}

source("KIRC_getRes.R")
kirc <- KIRC_getRes()
kirc2 <- .getBest(kirc)
kirc2 <- cbind(kirc2,"KIRC")

source("GBM_getRes.R")
gbm <- GBM_getRes()
gbm2 <- .getBest(gbm)
gbm2 <- cbind(gbm2,"GBM")

source("OV_getRes.R")
ov <- OV_getRes()
ov2 <- .getBest(ov)
ov2 <- cbind(ov2,"OV")

source("LUSC_getRes.R")
lusc <- LUSC_getRes()
lusc2 <- .getBest(lusc)
lusc2 <- cbind(lusc2,"LUSC")

comb <- rbind(kirc2,gbm2,lusc2,ov2)
comb <- as.data.frame(comb)
rownames(comb)<- NULL
comb[,1] <- as.numeric(as.character(comb[,1]))
colnames(comb) <- c("val","method","datatype","tumour")
dt <- format(Sys.Date(),"%y%m%d")
write.table(comb,file=sprintf("netDx_bestModel_%s.txt",dt),sep="\t",col=T,row=F,quote=F)

# convert into table form for comparison with Yuan et al.
x <- dcast(comb,tumour~datatype,value.var="val")
write.table(x,file=sprintf("netDx_perf_%s.txt",dt),sep="\t",col=T,row=F,quote=F)
