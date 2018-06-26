#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC/output/pathway80_170720

mkdir -p $localDir

# DELL
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathwaysOnly_top80_170719
outF=path80_170719.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

# DELL
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathwaysOnly_top80_2_170719
outF=path80_2_170719.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

mv *gz ${localDir}/.

tar xvfz KIRC_featSelPath_VM3_170719.tar.gz

