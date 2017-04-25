#!/bin/bash

#' count num net cols for all nets that go into making the GM database

inDir=/mnt/data2/BaderLab/PanCancer_GBM/input
netDir=/mnt/data2/BaderLab/PanCancer_GBM/output/ownTrain_170206/networks


echo "---------------------------------------"
echo "input data"
echo "---------------------------------------"
for f in ${inDir}/*core.txt;do
	echo $f
	head -n 1 $f | gawk '{print NF}';
done

echo "---------------------------------------"
echo "GM nets"
echo "---------------------------------------"
for f in ${netDir}/*.profile; do
	echo $f
	head -n 1 $f | gawk '{print NF}';
done


