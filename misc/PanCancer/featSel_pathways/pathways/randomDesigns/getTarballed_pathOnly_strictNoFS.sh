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
outF= VM1_sampledNets.txt
ssh $VM1 "cd $indir; ls randomPathNotFS_randomNets_*/SURVIVE*/networks/*.profile | xargs -L1 sh -c 'basename $1 .profile' dummy | sort -k1,1 > $outF"
scp $VM1:~/${outF} .

###### VM3
###indir=/home/spai/BaderLab/PanCancer_KIRC/output/
###outF=VM3_sampledNets.txt
###ssh $VM3 "cd $indir;
###
###scp $VM3:~/${outF} .
#####
###### VM5
###indir=/home/spai/BaderLab/PanCancer_KIRC/output/
###outF=KIRC_VM5randomStrictNoFS_170711.tar.gz
###ssh $VM5 "cd $indir; tar cvfz ~/${outF} randomPathNotFS_randomNets_*/predictions"
###scp $VM5:~/${outF} .
###
#### Dell
###indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output
###outF=KIRC_DellrandomStrictNoFS_170711.tar.gz
###ssh $DELL "cd $indir; tar cvfz ~/${outF} randomPathNotFS_randomNets_*/predictions"
###scp $DELL:~/${outF} .
###
###
#####
###mv *gz ${localDir}/.
###cd $localDir
###tar xvfz KIRC_VM1randomStrictNoFS_170711.tar.gz
###tar xvfz KIRC_VM3randomStrictNoFS_170711.tar.gz
###tar xvfz KIRC_VM5randomStrictNoFS_170711.tar.gz
###tar xvfz KIRC_DellrandomStrictNoFS_170711.tar.gz
