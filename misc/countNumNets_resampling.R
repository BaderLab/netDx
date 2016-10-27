#!/bin/bash

# count num nets in databases creating for 'pick best cutoff' phase

inDir=/home/spai/tmp/TCGA_BRCA_xprCNV/predictor

subtypes=(LumA other);
for k in {1..3};do
	for g in ${subtypes[@]};do
		ls ${inDir}/eval/part${k}/${g}/tmp/INTERACTIONS/*.txt | wc -l
	done
done
