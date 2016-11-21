# netDx output

This file describes netDx output files and formats

* [Overall directory tree](#overall_directory_tree)
* **File formats:** \[<a name="#profile">.profile</a>\][\_cont.txt][.NRANK][.PRANK][.query][pathway_CV_score.txt or pathway_cumTally.txt]

<a name="overall_directory_tree"></a>
## Overall directory tree
Assume the base directory for the dataset is `outRoot`.

### Example 1: BRCA LumA from gene expression and CNV: No resampling
Here we see the output in the case where the predictor is run once, with a single round of 10-fold cross validation. An example is the vignette building a predictor for Luminal A type of breast cancer from TCGA data.

See code snippet 1 below for the full set of output files generated in this example. Let `N` be the total number of patients in the dataset, and let `N'` be number of training samples.
GM = GeneMANIA, the algorithm used to integrate similarity networks and rank patients.
```
outRoot
	- profiles/  ### patient nets using training samples
		- *.profile : # input profiles to be converted to interaction nets
		- *_cont.txt: # interaction nets directly created by netDx
	- tmp/	# temporary files for GeneMANIA database (do not edit by hand)
		- GENES.txt: # list of patients in GM database
		- INTERACTIONS/*txt: # interaction nets as GM will see them
		- ... # (other files the user will not need)
	- dataset/ # location of indexed GeneMANIA database, not human-readable
	- LumA/ # results for class 1
		- GM_results 
			- CV_*.query # GM query file for a fold of CV 
			- CV_*.query-results.report.txt.PRANK # GM results file, patient ranking
			- CV_*.query-results.report.txt.NRANK # GM results file, network table
			- LumA_pathway_CV_score.txt # cumulative network score for 'LumA' predictive 
						    # value over entire CV
	- other/ # results for class 2
			- CV_*.query
			- CV_*.query-results.report.txt.PRANK
			- CV_*.query-results.report.txt.NRANK
			- other_pathway_CV_score.txt # cumulative network score for 'other' predictive
							# value over entire CV
```

## File formats

### cont.txt

### .NRANK

### .PRANK

<a name="#profile"></a>
### .profile

### .query
