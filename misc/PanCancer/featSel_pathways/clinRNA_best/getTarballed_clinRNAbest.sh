#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206
VM5=spai@192.168.81.122

localDir='/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/KIRC_clinRNA_best'

mkdir -p $localDir

# VM3
indir=/home/spai/BaderLab/PanCancer_KIRC/output/featSel_pathways_170519
outF1=KIRC_featSelPath_VM3_170502.tar.gz
ssh $VM3 "cd $indir; tar cvfz ~/${outF1} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM3:~/${outF1} .

# VM4
indir=/home/netdx/BaderLab/PanCancer_KIRC/output/featSel_pathways_170519
outF2=KIRC_featSelPath_VM4_170502.tar.gz
ssh $VM4 "cd $indir; tar cvfz ~/${outF2} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM4:~/${outF2} .

# VM5
indir=/home/spai/BaderLab/PanCancer_KIRC/output/featSel_pathways_170519
outF3=KIRC_featSelPath_VM5_170502.tar.gz
ssh $VM5 "cd $indir; tar cvfz ~/${outF3} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM5:~/${outF3} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz $outF1
tar xvfz $outF2
tar xvfz $outF3

