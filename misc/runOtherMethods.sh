#!/bin/bash
# run other methods to compare with netdx

# avg xpr by pathway
#Rscript enet_avgXprByPathway.R ElasticNet
#Rscript enet_avgXprByPathway.R RandomForest
# genes in pathways
Rscript enet_genesInPathways.R ElasticNet
Rscript enet_genesInPathways.R RandomForest


