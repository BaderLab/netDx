#!/bin/bash
# run on mac, collects nets used in each round of random sampling

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206
VM5=spai@192.168.81.122
DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/KIRC_pathOnly_randomStrictNoFS_170711

mkdir -p $localDir
##
### VM1
indir=/mnt/data2/BaderLab/PanCancer_KIRC/output/
outF=VM1_sampledNets.txt
ssh $VM1 "cd $indir; ls randomPathNotFS_randomNets_*/SURVIVE*/networks/*.profile  > ~/$outF"
scp $VM1:~/${outF} .

###### VM3
indir=/home/spai/BaderLab/PanCancer_KIRC/output/
outF=VM3_sampledNets.txt
ssh $VM3 "cd $indir; ls randomPathNotFS_randomNets_*/SURVIVE*/networks/*.profile  > ~/$outF"
scp $VM3:~/${outF} .

#### Dell
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output
outF=Dell_sampledNets.txt
ssh $DELL "cd $indir; ls randomPathNotFS_randomNets_*/SURVIVE*/networks/*.profile  > ~/$outF"
scp $DELL:~/${outF} .

