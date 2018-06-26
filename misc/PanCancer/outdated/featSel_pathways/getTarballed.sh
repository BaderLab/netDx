#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206

localDir=/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC/output/KIRC_featSel_pathway_170426

mkdir -p $localDir

# VM1
indir=/mnt/data2/BaderLab/PanCancer_KIRC/output/featSel_pathways_170426
outF=KIRC_featSelPath_VM1_170426.tar.gz
ssh $VM1 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM1:~/${outF} .

# VM3
indir=/home/spai/BaderLab/PanCancer_KIRC/output/featSel_pathways_170426
outF=KIRC_featSelPath_VM3_170426.tar.gz
ssh $VM3 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM3:~/${outF} .

# VM4
indir=/home/netdx/BaderLab/PanCancer_KIRC/output/featSel_pathways_VM4_170426
outF=KIRC_featSelPath_VM4_170426.tar.gz
ssh $VM4 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM4:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_featSelPath_VM1_170426.tar.gz
tar xvfz KIRC_featSelPath_VM3_170426.tar.gz
tar xvfz KIRC_featSelPath_VM4_170426.tar.gz

