This file describes the content of detailed output by the netDx algorithm.

Assume the base directory for the dataset is `outRoot`.

## Example 1: BRCA LumA from gene expression and CNV: No resampling
Here we see the output in the case where the predictor is run once, with a single round of 10-fold cross validation. An example is the vignette building a predictor for Luminal A type of breast cancer from TCGA data.

See code snippet 1 below for the full set of output files generated in this example. 


```
outRoot
	- LumA/
		- GM_results
			- CV_*.query
			- CV_*.query-results.report.txt.PRANK
			- CV_*.query-results.report.txt.NRANK
			- LumA_pathway_CV_score.txt
	- other/ 
			- CV_*.query
			- CV_*.query-results.report.txt.PRANK
			- CV_*.query-results.report.txt.NRANK
			- other_pathway_CV_score.txt
	- profiles/
	- tmp/
	- dataset/
```
In this example, we:
1. split data into train and test
### Patient nets from all training samples
`outRoot/tmp`
	* Build a master GM database with these input nets (`outRoot/tmp` and `outRoot/dataset`
3. apply feature selection once for each of the two classes, `LumA` and `other` (results in the `LumA` and `other` directory respectively)
4. following feature selection, we classify test patients by ranking them against two class-specific databa

Then after running netDx, the data directory looks like this:


