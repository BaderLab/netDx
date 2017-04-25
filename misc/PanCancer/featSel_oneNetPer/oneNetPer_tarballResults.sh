#!/bin/bash
# tarballs all prediction results and pathway scores.

# OV - VM2
cd /home/ahmad/tcga_datasets/OV/output/featSel_oneNetPer_170425 
outDir=/home/spai/BaderLab/PanCancer_common
tar cvfz ${outDir}/OV_oneNetPer_170425.tar.gz rng*/*/predictionResults.txt rng*/*/SURVIVE*/GM_results/*pathway_CV_score.txt

# KIRC - VM1

