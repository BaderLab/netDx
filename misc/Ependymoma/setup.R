#' 

inDir <- "~/Dropbox/netDx/BaderLab/2017_Ependymoma/input"
inFile <- sprintf("%s/GSE27279.Rdata",inDir)
sampFile <- sprintf("%s/WittPfister_TableS2.txt",inDir)

load(inFile)
samp <- read.delim(sampFile,sep="\t",h=T,as.is=T)
