#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206
VM5=spai@192.168.81.122

#### get KIRC
###localDir=/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC/output/KIRC_oneNetPer_normDiff_170518
###mkdir -p $localDir
###dir1=/mnt/data2/BaderLab/PanCancer_KIRC/output/oneNetPer_normDiff_170517
###dir2=/mnt/data2/BaderLab/PanCancer_KIRC/output/oneNetPer_normDiff_170518
###outF1=KIRC_normDiff_170518_part1.tar.gz
###ssh $VM1 "cd $dir1; tar cvfz ~/${outF1} rng*/clin*/predictionResults.txt rng*/clin*/SURVIVE*/*/*CV_score.txt rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt"  
###outF2=KIRC_normDiff_170518_part2.tar.gz
###ssh $VM1 "cd $dir2; tar cvfz ~/${outF2} rng*/clin*/predictionResults.txt rng*/clin*/SURVIVE*/*/*CV_score.txt rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt"  
###
###scp $VM1:~/${outF1} .
###scp $VM1:~/${outF2} .
###
###mv *gz ${localDir}/.
###cd $localDir
###tar xvfz $outF1
###tar xvfz $outF2
###

#### get GBM
localDir=/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_GBM/output/GBM_oneNetPer_normDiff_170518
mkdir -p $localDir
dir1=/home/spai/BaderLab/PanCancer_GBM/output/oneNetPer_normDiff_170517
dir2=/home/spai/BaderLab/PanCancer_GBM/output/oneNetPer_normDiff_170518
outF1=GBM_normDiff_170518_part1.tar.gz
ssh $VM3 "cd $dir1; tar cvfz ~/${outF1} rng*/clin*/predictionResults.txt rng*/clin*/SURVIVE*/*/*CV_score.txt rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt"  
outF2=GBM_normDiff_170518_part2.tar.gz
ssh $VM3 "cd $dir2; tar cvfz ~/${outF2} rng*/clin*/predictionResults.txt rng*/clin*/SURVIVE*/*/*CV_score.txt rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt"  

scp $VM3:~/${outF1} .
scp $VM3:~/${outF2} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz $outF1
tar xvfz $outF2

