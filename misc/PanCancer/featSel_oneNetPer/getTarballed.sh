#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206

localDir=/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC/output/KIRC_oneNetPer_170425

mkdir -p $localDir


# VM4 - gets 89-100
indir=/home/netdx/BaderLab/PanCancer_KIRC/output/featSel_oneNetPer_170426
outF=KIRC_oneNetPer_VM4_170426.tar.gz
ssh $VM4 "cd $indir; tar cvfz ~/${outF} rng*/*/predictionResults.txt rng*/*/SURVIVE*/*/*CV_score.txt" 
scp $VM4:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_oneNetPer_VM4_170426.tar.gz

