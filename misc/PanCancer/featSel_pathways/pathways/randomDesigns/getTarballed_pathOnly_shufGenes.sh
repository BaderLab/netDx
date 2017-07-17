#!/bin/bash
# run on mac, collects nets used in each round of random sampling

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206
VM5=spai@192.168.81.122
DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/KIRC_pathOnly_randomShuf_170712

mkdir -p $localDir
##
### VM1
indir=/mnt/data2/BaderLab/PanCancer_KIRC/output/
outF=KIRC_VM1randomShuf_070712.tar.gz
ssh $VM1 "cd $indir; tar cvfz ~/${outF} randomPathShuf_randomNets_*/predictions"
scp $VM1:~/${outF} .

######### VM3
indir=/home/spai/BaderLab/PanCancer_KIRC/output/
outF=KIRC_VM3randomShuf_070712.tar.gz
ssh $VM3 "cd $indir; tar cvfz ~/${outF} randomPathShuf_randomNets_*/predictions"
scp $VM3:~/${outF} .

# Dell
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output
outF=KIRC_DellrandomShuf_070712.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} randomPathShuf_randomNets_*/predictions"
scp $DELL:~/${outF} .
###
mv *gz ${localDir}/.
cd $localDir
tar xvfz KIRC_VM1randomShuf_070712.tar.gz
tar xvfz KIRC_VM3randomShuf_070712.tar.gz
tar xvfz KIRC_DellrandomShuf_070712.tar.gz
