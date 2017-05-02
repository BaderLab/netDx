#!/bin/bash
# run on mac, remotely tarballs results and downloads
# random results for KIRC

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206

localDir=/Users/shraddhapai/Documents/Research/BaderLab/2017_TCGA_KIRC/output/KIRC_oneClinNet_pathway_random_170501

mkdir -p $localDir

# VM1
indir=/mnt/data2/BaderLab/PanCancer_KIRC/output/
outF=KIRC_featSelPath_VM1_170501.tar.gz
ssh $VM1 "cd $indir; tar cvfz ~/${outF} randomNets_forceClin*/predictions"
scp $VM1:~/${outF} .

# VM4
indir=/home/netdx/BaderLab/PanCancer_KIRC/output/
outF=KIRC_featSelPath_VM4_170501.tar.gz
ssh $VM4 "cd $indir; tar cvfz ~/${outF} randomNets_forceClin*/predictions"
scp $VM4:~/${outF} .

mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_featSelPath_VM1_170501.tar.gz
tar xvfz KIRC_featSelPath_VM4_170501.tar.gz

