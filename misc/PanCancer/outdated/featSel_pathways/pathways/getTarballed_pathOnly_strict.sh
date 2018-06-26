#!/bin/bash
# run on mac, remotely tarballs results and downloads
# random results for KIRC

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206
VM5=spai@192.168.81.122
DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_TCGA_KIRC/output/KIRC_pathOnly_randomStrict_170620

mkdir -p $localDir
##
### VM1
##indir=/mnt/data2/BaderLab/PanCancer_KIRC/output/
##outF=KIRC_VM1randomStrict_170620.tar.gz
##ssh $VM1 "cd $indir; tar cvfz ~/${outF} pathOnly_randomNets*170607/predictions"
##scp $VM1:~/${outF} .
##
### VM3
###indir=/home/spai/BaderLab/PanCancer_KIRC/output/
###outF=KIRC_VM3randomStrict_170620.tar.gz
###ssh $VM3 "cd $indir; tar cvfz ~/${outF} pathOnly_randomNets*170607/predictions"
###scp $VM3:~/${outF} .
##
### VM5
##indir=/home/spai/BaderLab/PanCancer_KIRC/output/
##outF=KIRC_VM5randomStrict_170620.tar.gz
##ssh $VM5 "cd $indir; tar cvfz ~/${outF} pathOnly_randomNets*170607/predictions"
##scp $VM5:~/${outF} .

# VM5 - part 2
indir=/home/spai/BaderLab/PanCancer_KIRC/output/
outF=KIRC_VM5randomStrict_170709.tar.gz
ssh $VM5 "cd $indir; tar cvfz ~/${outF} pathOnly_randomNets*170707/predictions"
scp $VM5:~/${outF} .

# Dell
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output
outF=KIRC_DellrandomStrict_170709.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} pathOnly_randomNets*170707/predictions"
scp $DELL:~/${outF} .


##
mv *gz ${localDir}/.
cd $localDir
##tar xvfz KIRC_VM1randomStrict_170620.tar.gz
##tar xvfz KIRC_VM3randomStrict_170620.tar.gz
##tar xvfz KIRC_VM5randomStrict_170620.tar.gz
tar xvfz KIRC_VM5randomStrict_170709.tar.gz
tar xvfz KIRC_DellrandomStrict_170709.tar.gz
