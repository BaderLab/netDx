#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206
VM5=spai@192.168.81.122

localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output/pathway_170502_forNetScoreTest

mkdir -p $localDir

# VM1
indir=/mnt/data2/BaderLab/PanCancer_KIRC/output/pathwaysOnly_170502
outF=KIRC_featSelPath_VM1_170502.tar.gz
ssh $VM1 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/GM_results rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM1:~/${outF} .

# VM3
indir=/home/spai/BaderLab/PanCancer_KIRC/output/pathwaysOnly_VM3_170502
outF=KIRC_featSelPath_VM3_170502.tar.gz
ssh $VM3 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/GM_results rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM3:~/${outF} .

# VM4
indir=/home/netdx/BaderLab/PanCancer_KIRC/output/pathwaysOnly_170502
outF=KIRC_featSelPath_VM4_170502.tar.gz
ssh $VM4 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/GM_results rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM4:~/${outF} .

# VM5
indir=/home/spai/BaderLab/PanCancer_KIRC/output/pathwaysOnly_VM5_v2_170502
outF=KIRC_featSelPath_VM5_170502.tar.gz
ssh $VM5 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/GM_results rng*/SURVIVE*/*/*CV_score.txt" 
scp $VM5:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_featSelPath_VM1_170502.tar.gz
tar xvfz KIRC_featSelPath_VM3_170502.tar.gz
tar xvfz KIRC_featSelPath_VM4_170502.tar.gz
tar xvfz KIRC_featSelPath_VM5_170502.tar.gz

