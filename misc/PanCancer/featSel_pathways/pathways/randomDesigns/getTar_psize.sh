#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output/pathSize_170808
mkdir -p $localDir

indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/changePathSize7_170808
outF=KIRC_changePathSize_170808.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} pSize*/numG*/rng*/predictionResults.txt"
scp $DELL:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_changePathSize_170808.tar.gz

