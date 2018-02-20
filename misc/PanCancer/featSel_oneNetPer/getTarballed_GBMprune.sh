#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

DELL=shraddhapai@192.168.81.205

localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_GBM/output/pruned_180204

mkdir -p $localDir

indir=/home/shraddhapai/BaderLab/2017_PanCancer/GBM/output/prune_180204
outF=GBM_pruned_180204.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/*/cutoff9/predictionResults.txt rng*/*/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz $outF

