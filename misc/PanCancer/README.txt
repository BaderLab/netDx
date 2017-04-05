This folder contains code to run the PanCancer survival related predictor
tests.

## netDx predictor runs - one net per datatype
* OV_classify.R ; OV_classify_ownTrain.R
* LUSC_classify.R ; LUSC_classify_ownTrain.R

## netDx predictor runs - feature selection
* GBM_featSel.R ; GBM_perf.R
* LUSC_featSel.R ; LUSC_perf.R
	* LUSC_featSel_clinRNA.R ; LUSC_featSel_clinRPPA.R 
	* LUSC_featSel_somMut.R
* OV_featSel_somMut.R ; OV_perf.R
* Get mean+SEM AUC for all cancer types: PanCancer_featSel_plotPerf.R

## Making enrichment maps
* Base: PanCancer_featSel_writeConsensusMap.R
	* PanCancer_callConsMapMaker.R

## Integrated PSN
* GBM_getPSN.R

