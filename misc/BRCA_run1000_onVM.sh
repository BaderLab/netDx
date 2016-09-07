#!/bin/bash

# writes and qsubs many jobs for BRCA prediction
resDir=/home/spai/tmp/TCGA_BRCA_runMany

numRuns=250
for (( k=120; k<=$numRuns; k++ )); do
    outDir=${resDir}/R${k}
	chmod -R u+w $outDir; rm -r $outDir;
	mkdir $outDir
    seed1=$((10 * $k))
    seed2=$((15 * $k))
    echo "${k}: Seeds are ${seed1} and ${seed2}"

    Rscript BRCA_run1000_scinet.R $outDir $seed1 $seed2
    #chmod -R u-w $outDir
done
