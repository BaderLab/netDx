#!/bin/bash

# writes and qsubs 1000 jobs for BRCA prediction
jobDir=/scratch/g/gbader/spai/BaderLab/BRCA/jobs
resDir=/scratch/g/gbader/spai/BaderLab/BRCA/results
scriptdir=/home/g/gbader/spai/software/netDx/misc
wtime="2:00:00"

numRuns=1

for (( k=1; k<=$numRuns; k++ )); do
    jobFile=${jobDir}/job${k}.sh
    outDir=${resDir}/R${k}
    seed1=$((10 * $k))
    seed2=$((15 * $k))
    echo "${k}: Seeds are ${seed1} and ${seed2}"

    cat /dev/null > $jobFile;
    echo "#!/bin/bash" >> $jobFile
    echo "#PBS -l nodes=1:ppn=8,walltime=${wtime}" >> $jobFile
    echo "#PBS -N BRCA.Run${k}" >> $jobFile
    echo "cd $scriptdir;"           >> $jobFile
    echo ""                         >> $jobFile
    echo "module load java/6.0 gcc intel/15.0.6 R/3.2.3" >> $jobFile
    echo "Rscript BRCA_run1000_scinet.R $outDir $seed1 $seed2" >> $jobFile
    echo "chmod -R u+w $outDir" >> $jobFile # protect results!

    chmod u+x $jobFile
	curd=`pwd`
	cd $jobDir
	qsub $jobFile
	cd $curd
done
