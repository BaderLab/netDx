#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output/simRank_170719/designD_rmFSgenes_170717
mkdir -p $localDir

indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathRandom_designD_rmFSgenes_170717
outF=KIRC_randomD_noFS_170717.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/SURVIVE*/*PRANK"
scp $DELL:~/${outF} .

indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathRandom_designD_rmFSgenes2_170717
outF=KIRC_randomD_noFS2_170717.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/SURVIVE*/*PRANK"
scp $DELL:~/${outF} .


mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_randomD_noFS2_170717.tar.gz
tar xvfz KIRC_randomD_noFS_170717.tar.gz

