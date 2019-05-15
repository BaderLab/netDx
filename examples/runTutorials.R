
runTut <- function(tFile) {
	logFile <- sprintf("%s.log",tFile)
	sink(logFile,split=TRUE)
	print(Sys.time())
tryCatch({
	
	rmarkdown::render(tFile)
},error=function(ex) {
	print(ex)
},finally={
	cat(sprintf("%s done\n", basename(tFile)))
	print(Sys.time())
	sink(NULL)	
})
}

tutFiles <- c("Medulloblastoma.Rmd","BreastCancer.Rmd","BuildPredictor.Rmd")
for (tFile in tutFiles) runTut(tFile)

