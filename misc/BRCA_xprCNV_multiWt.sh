#!/bin/bash

# run netDx predictor for BRCA xpr+cnv once on hungabee.westgrid.ca
### Loop over different filter_WtSum values

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

#module load application/jdk/1.8.0_92


# --------------

outRoot=/home/spai/tmp/BRCA_xprCNV_multiSeed
mkdir -p $outRoot

wt=(100 90 75 50);

for k in {200..204}; do
	seed1=$k
	seed2=$(( k * 10 ))
	echo "seed1=${k} ; seed2=${k}"
	omain=${outRoot}/R${k}
	echo $omain
	mkdir -p $omain
	## loop over different net weight inclusions
	for w in ${wt[@]};do
		odir=${omain}/wt${w}
		mkdir -p $odir
		Rscript BRCA_xprCNV_resample_changeWtSum.R \
			$odir $seed1 $seed2 $w
		rm -r ${odir}/predictor
		chmod u-w ${odir}/*.*
	done
done


