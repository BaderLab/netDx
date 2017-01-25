#' feature selection for DREAM AD Question 3
#' here we start with the input files processed by GuanLab, one of the
#' winning teams of DREAM 2014 AD challenge.
#' 
#' we perform feature selection over the three diagnostic groups and 
#' then do a final classification.

numCores <- 8L
GMmemory <- 4L
trainProp <- 0.67

require(netDx)
inFile <- "/home/spai/BaderLab/DREAM_AD/input/GuanLab/C3/ADNI_Training_Q3_new.csv_matlab.csv"
outRoot <- "/home/spai/BaderLab/DREAM_AD/output"

dt <- format(Sys.Date(),"%y%m%d")
outDir <- sprintf("%s/featSel_%s",dt)

# ----------------------------------------------------------------

dat <- read.delim(inFile,sep=",",h=F,as.is=T)

# imaging networks have correlation by similarity
# clinical/genetic networks have similarity by normalized distance.

