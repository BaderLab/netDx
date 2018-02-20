#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

DELL=shraddhapai@192.168.81.205

localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output/pruned_180204

mkdir -p $localDir

indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pruned_180204
outF=KIRC_pruned_180204.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/*/predictionResults.txt rng*/*/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz $outF

