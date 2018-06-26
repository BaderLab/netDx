#!/bin/bash
# run on mac, remotely tarballs results and downloads
# random results for KIRC

DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/pathways_AltClass_noFS_170721

mkdir -p $localDir

# Dell
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathways_AltClass_noFS_170721
outF=KIRC_altClass_noFS_170721.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

cd $localDir
tar xvfz KIRC_altClass_noFS_170721.tar.gz
