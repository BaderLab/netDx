#!/bin/bash

sig_val=(0.05 0.1 0.3 1);
for cursig in ${sig_val[@]}; do
	echo $cursig
	screen -d -m bash -c "Rscript GBM_univar_rbf.R rbf $cursig"
	sleep 5
done

