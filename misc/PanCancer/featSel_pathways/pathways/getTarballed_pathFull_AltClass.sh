#!/bin/bash
# run on mac, remotely tarballs results and downloads
# random results for KIRC

DELL=shraddhapai@192.168.81.215
VM1=spai@192.168.81.125
VM3=spai@192.168.81.249

localDir=/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/pathFull_AltClass_170721

mkdir -p $localDir

#### Dell
###indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathwayFull_AltClass_170721
###outF=KIRC_altClass_170721.tar.gz
###ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
###scp $DELL:~/${outF} .
###
#### Dell - 2
###indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/pathwayFull_AltClass_2_170721
###outF=KIRC_altClass2_170721.tar.gz
###ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
###scp $DELL:~/${outF} .
###
#### VM1
###indir=/mnt/data2/BaderLab/PanCancer_KIRC/output/pathwayFull_AltClass_170721
###outF=KIRC_altClass_VM1_170721.tar.gz
###ssh $VM1 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
###scp $VM1:~/${outF} .
###
#### VM3
###indir=/home/spai/BaderLab/PanCancer_KIRC/output/pathwayFull_AltClass_170721
###outF=KIRC_altClass_VM3_170721.tar.gz
###ssh $VM3 "cd $indir; tar cvfz ~/${outF} rng*/predictionResults.txt rng*/SURVIVE*/*/*CV_score.txt" 
###scp $VM3:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_altClass_170721.tar.gz
tar xvfz KIRC_altClass2_170721.tar.gz
tar xvfz KIRC_altClass_VM1_170721.tar.gz
tar xvfz KIRC_altClass_VM3_170721.tar.gz
