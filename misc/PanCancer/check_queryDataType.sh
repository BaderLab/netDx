#!/bin/bash

# check that query is indeed using the datatypes we think it is

inRoot=/mnt/data2/BaderLab/PanCancer_GBM/output/ownTrain_170206/run45

for f in ${inRoot}/*/SURVIVEYES_testQuery;do
	echo $f
	sed -n '3p' $f
done

