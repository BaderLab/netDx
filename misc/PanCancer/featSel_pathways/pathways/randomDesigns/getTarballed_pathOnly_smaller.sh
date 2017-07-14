#!/bin/bash
# run on mac, collects nets used in each round of random sampling

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
DELL=shraddhapai@192.168.81.215

#### SMALLER

#### dell - smaller
###indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output
###outF=KIRC_DellrandomSmaller_170713.tar.gz
###ssh $DELL "cd $indir; tar cvfz ~/${outF} randomPathSmaller_*/predictions"
###scp $DELL:~/${outF} .
###localDir=/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/KIRC_pathOnly_randomSmaller_170713
###mkdir -p $localDir
###mv *gz ${localDir}/.
###cd $localDir
###tar xvfz $outF

# dell - smallest
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output
outF=KIRC_DellrandomSmallest_170713.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} randomPathSmallest_*/predictions"
scp $DELL:~/${outF} .
localDir=/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/KIRC_pathOnly_randomSmallest_170713
mkdir -p $localDir
mv *gz ${localDir}/.
cd $localDir
tar xvfz $outF

