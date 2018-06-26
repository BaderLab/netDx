#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output/randomD_170713
mkdir -p $localDir

indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathRandom_designD_170713
outF=KIRC_randomD_170713.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_featSelPath_VM1_170502.tar.gz

