#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

DELL=shraddhapai@192.168.81.215

#### DELL
###localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output/pathway90_170720
###mkdir -p $localDir
###indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathwaysOnly_topX_170720/top90
###outF=path90_170720.tar.gz
###ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
###scp $DELL:~/${outF} .
###mv *gz ${localDir}/.
###cd $localDir
###tar xvfz path90_170720.tar.gz

# DELL
localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output/pathway95_170720
mkdir -p $localDir
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathwaysOnly_topX_2_170720/top95
outF=path95_2_170720.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .
mv *gz ${localDir}/.

cd $localDir
tar xvfz path95_2_170720.tar.gz
