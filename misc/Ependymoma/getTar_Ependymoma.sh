#!/bin/bash

srcDir=/home/shraddhapai/BaderLab/2017_Ependymoma/output/Epen_180125/pred
tgtDir=/Users/shraddhapai/Dropbox/netDx/BaderLab/2017_Ependymoma/output/Epen_180125/pred
ssh $DELL "cd $srcDir; tar cvfz res.tar.gz rng*/predictionResults.txt rng*/*/GM_results/*pathway_CV_score.txt inputNets.txt log.txt rng1/*/tmp rng2/*/tmp"
mkdir -p $tgtDir
scp $DELL:${srcDir}/res.tar.gz $tgtDir
cd $tgtDir
tar xvfz res.tar.gz

