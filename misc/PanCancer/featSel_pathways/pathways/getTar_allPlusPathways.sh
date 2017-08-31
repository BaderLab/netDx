#!/bin/bash
# run on mac, remotely tarballs results and downloads
# KIRC pathway predictor

VM1=spai@192.168.81.125
VM3=spai@192.168.81.249
VM4=netdx@192.168.81.206
VM5=spai@192.168.81.122
DELL=shraddhapai@192.168.81.215

localDir=/Users/shraddhapai/DropBox/netDx/BaderLab/2017_TCGA_KIRC/output/AllPlusPathways_170811
mkdir $localDir

# DELL
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/AllPlusPathways_170811
outF=AllPlusPathways_Dell_170811.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

# DELL2
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/AllPlusPathways2_170811
outF=AllPlusPathways_Dell2_170811.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

# DELL3
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/AllPlusPathways3_170811
outF=AllPlusPathways_Dell3_170811.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

# DELL4
indir=/home/shraddhapai/BaderLab/PanCancer_KIRC/output/AllPlusPathways4_170813
outF=AllPlusPathways_Dell4_170813.tar.gz
ssh $DELL "cd $indir; tar cvfz ~/${outF} rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt" 
scp $DELL:~/${outF} .

# VM1
indir=/mnt/data2/BaderLab/PanCancer_KIRC/output/AllPlusPathways2_170811
outF=AllPlusPathways_VM1_170811.tar.gz
ssh $VM1 "cd $indir; tar cvfz ~/${outF} rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt" 
scp $VM1:~/${outF} .

# VM3
indir=/home/spai/BaderLab/PanCancer_KIRC/output/AllPlusPathways2_170811
outF=AllPlusPathways_VM3_170811.tar.gz
ssh $VM3 "cd $indir; tar cvfz ~/${outF} rng*/all/predictionResults.txt rng*/all/SURVIVE*/*/*CV_score.txt" 
scp $VM3:~/${outF} .


mv *gz ${localDir}/.
cd $localDir
tar xvfz AllPlusPathways_VM1_170811.tar.gz
tar xvfz AllPlusPathways_VM3_170811.tar.gz
#tar xvfz AllPlusPathways_VM4_170811.tar.gz
tar xvfz AllPlusPathways_Dell4_170813.tar.gz
tar xvfz AllPlusPathways_Dell3_170811.tar.gz
tar xvfz AllPlusPathways_Dell2_170811.tar.gz
tar xvfz AllPlusPathways_Dell_170811.tar.gz
