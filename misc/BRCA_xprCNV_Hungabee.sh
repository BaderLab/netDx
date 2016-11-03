#!/bin/bash

# run netDx predictor for BRCA xpr+cnv once on hungabee.westgrid.ca

### from https://www.westgrid.ca/support/quickstart/hungabee:
### It is usually easiest to keep pmem=8190mb and adjust procs to insure that there is sufficient total memory available. For example, if you need 2000000mb of memory, 2000000/8192=244.1. However, procs needs to be a multiple of 16 so procs should be 256.

# SP's note: the job below is asking for 50Gb RAM, i.e. procs is set to
# (50*1024)Mb / (8190mb) per proc = 6.25 procs. Round up to 16.

#PBS -l procs=16
#PBS -l pmem=8190mb
#PBS -l walltime=2:30:00
#PBS -N BRCA_test
#PBS -m abe
#PBS -M shraddha.pai@utoronto.ca

module load application/jdk/1.8.0_92

Rscript=/home/spai/software/R-3.2.5/bin/Rscript
cd /home/spai/software/netDx/misc

# --------------

outRoot=/data/spai/BaderLab/netDx/BRCA_xprCNV

for k in 1; do
	seed1=$k
	seed2=$(( k * 10 ))
	echo "seed1=${k} ; seed2=${k}"
	odir=${outRoot}/R${k}
	mkdir $odir
	$Rscript BRCA_xprCNV_resample.R $odir $seed1 $seed2
done


