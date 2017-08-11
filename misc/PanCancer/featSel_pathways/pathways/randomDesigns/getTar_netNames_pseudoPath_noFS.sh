#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_PanCancer_Survival/randomD_pseudoPath_noPathGenes_170804
mkdir $localDir

indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/randomD_pseudoPath_noPathGenes2_170804
outF=pseudoNets2_170804.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/SURVIVE*/netNames.txt"
scp $DELL:~/${outF} .

indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/randomD_pseudoPath_noPathGenes_170804
outF=pseudoNets_170804.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/SURVIVE*/netNames.txt"
scp $DELL:~/${outF} .


mv *gz ${localDir}/.
cd $localDir
tar xvfz pseudoNets2_170804.tar.gz
tar xvfz pseudoNets_170804.tar.gz

