# netDx output

This file describes netDx output files and formats

* [Overall directory tree](#overall_dir)
	* [Example 1: Luminal A prediction from gene expresion and CNV: No resampling](#ex1)
	* [Example 2: Predictor with three resamplings](#ex2)
	* [Example 3: Predictor with nested cross-validation](#ex3)
* **File formats:** \[<a href="#intfile">\*\_cont.txt</a>][<a href="#nrank">\*.NRANK</a>\]\[<a href="#prank">\*.PRANK</a>\]\[<a href="#cvscore">pathway_CV_score.txt or pathway_cumTally.txt</a>\]\[<a href="#profile">\*.profile</a>\]\[<a href="#query">\*.query</a>\]

<a name="overall_dir"></a>
## Overall directory tree
Assume the base directory for the dataset is `outRoot`.

<a name="ex1"></a>
### Example 1: Luminal A predictfrom gene expression and CNV: No resampling
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

<a name="ex2"></a>
### Example 2: Predictor with 3 resamplings
The directory structure is similar to the previous example, with two notable differences:
1. results for each resampling are stored in their own subdirectory (e.g. `part1/`, `part2/`, etc.,)
2. there are different main directories for training vs test data.  Training results are stored in `feature_select/` while test are stored in `networks_test/`.

In this example, we have a binary classifier and feature selection is only run for the class in question (rather than for A and not-A).

```
outRoot
	- feature_select/ # training results
		- networks_orig/     # similarity nets directly made by netDx
			- *_cont.txt # input similarity nets
		- part1/  # results for first round of resampling
			- GM_results/ # GeneMANIA input/output files
				- CV_*.query  # query file for a fold of cross-validation
				- CV_*.query-results.report.txt.PRANK # GM output: patient ranking
				- CV_*.query-results.report.txt.NRANK # GM output: network weights
				- pathwayScore.txt # cumulative net score for current resampling
			- networks/   # interaction nets
			- networksCliqueFilt/ # subset of networks that pass clique filtering
			- cliqueFilterNets.stats.txt # net score for clique-filtering
			- tmp/ # input files to create GM database
			- dataset/ # indexed GM database on which queries are run (not human-readable)
		- part2/ # same structure as part1
		- part3/ # same structure as part1
		- pathway_cumTally.txt # cumulative net score for predictive class	
```
<a name="ex3"></a>
### Example 3: Predictor with nested cross-validation
In this design, the predictor has two loops for feature selection
#. An outer loop where the data are split into several (say K1) train/blind test groups
#. An inner loop where 10-fold cross validation is performed on each split.
The result is that there are K1 sets of network scores, each ranging from 0 to 10. There are also K1 sets of blind test predictions from which the mean and variation of predictor performance can be measured.

The output format recommended by netDx, and expected by downstream analysis code is as follows:
```
dataset_yymmdd/
  + rng1/
    + tmp/       # directory created by netDx, containing input data for GeneMANIA database
    + networks/  # PSN created by calls to makePSN_NamedMatrix()
    +-- Class1
       + tmp/
       + networks/                               # networks for test classification for this split
       + GM_results/                             # results of inner loop (10-fold CV)
       + Class1_pathway_CV_score.txt             # network scores for inner CV fold
       +--- CV_1.query                           # query for CV fold
       +--- CV_1.query-results.report.txt.NRANK  # network weights for CV fold
       +--- CV_1.query-results.report.txt.PRANK  # patient sim ranking for CV fold (not used)
       ...
       +--- CV_10.query                           # query for CV fold
       +--- CV_10.query-results.report.txt.NRANK  # network weights for CV fold
       +--- CV_10.query-results.report.txt.PRANK  # patient sim ranking for CV fold (not used)
    +-- Class2
    + predictionResults.txt  # test predictions for this train/test split
  + rng2/
  ...
  + rngK1/
```

## File formats
<a name="intfile"></a>
#### cont.txt
An individual patient similarity network, one per input feature.
For example, if features are pathways, each net is named after a pathway (e.g. WNT_SIGNALING_NETWORK_cont.txt, XENOBIOTICS_cont.txt).
Tab-delimited file without a header row; contains three columns: source patient, target patient, similarity weight.

Here is an example of TCGA patients that share CNVs in the same pathway. Patients without co-occurring CNVs are excluded.
```
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A2.A0CL.01A.11R.A115.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A2.A0D1.01A.11R.A034.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A2.A0YF.01A.21R.A109.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A7.A0D9.01A.31R.A056.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A8.A07O.01A.11R.A00Z.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A8.A081.01A.11R.A00Z.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A8.A083.01A.21R.A00Z.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A8.A08A.01A.11R.A00Z.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.A8.A08H.01A.21R.A00Z.07    1
TCGA.A1.A0SO.01A.22R.A084.07    TCGA.AN.A0AL.01A.11R.A00Z.07    1
```
Here is an example of a continuously-weighted network (e.g. based on gene expression correlation):
```
10      2       0.20394
15      2       0.23392
16      2       0.26524
27      2       0.2171
28      2       0.21943
54      2       0.18568
63      2       0.18642
74      2       0.15964
76      2       0.16618
78      2       0.17062
```

<a name="nrank"></a>
#### .NRANK
The part of a GeneMANIA query output that contains predictive weights of input networks.
Tab-delimited file with a header row.
The main columns of interest are 'Network' and 'Weight'.
```
Network Group   Network Weight  Title   Authors Year    Publication     PMID    URL     Processing Method       Interactions    Source  Source URL     Tags
dummy_group             100.00
        RHO_GTPASES_ACTIVATE_CIT.profile        3.46                                                            0
        NICOTINE_PHARMACODYNAMICS_PATHWAY.profile       3.02                                                            0
        BIOCARTA_RANMS_PATHWAY.profile  2.80                                                            0
        FOXM1_TRANSCRIPTION_FACTOR_NETWORK.profile      2.63                                                            0
        KINESINS.profile        2.56                                                            0
        ASSOCIATION_OF_LICENSING_FACTORS_WITH_THE_PRE-REPLICATIVE_COMPLEX.profile       2.44                                                           0

        DE_NOVO_PYRIMIDINE_DEOXYRIBONUCLEOTIDE_BIOSYNTHESIS.profile     2.36                                                            0
        AURORA_B_SIGNALING.profile      2.17 
```

<a name="prank"></a>
#### .PRANK
The part of a GeneMANIA query output that contains a table of patient rankings.
Tab-delimited file with a header row. Contains three columns:

1. Gene: patient ID
2. Score: blank for patients submitted as the query. In decreasing order for patients. Higher scoring patients are more similar to the query.
3. Description: (ignore)

```
$  head CV_3.query-results.report.txt.PRANK 
Gene    Score   Description
TCGA.A1.A0SH.01A.11R.A084.07            TCGA.A1.A0SH.01A.11R.A084.07
TCGA.B6.A0RN.01A.12R.A084.07            TCGA.B6.A0RN.01A.12R.A084.07
TCGA.BH.A0DP.01A.21R.A056.07            TCGA.BH.A0DP.01A.21R.A056.07
TCGA.AN.A0FS.01A.11R.A034.07            TCGA.AN.A0FS.01A.11R.A034.07
TCGA.E2.A1B4.01A.11R.A12P.07            TCGA.E2.A1B4.01A.11R.A12P.07
TCGA.A2.A0EO.01A.11R.A034.07            TCGA.A2.A0EO.01A.11R.A034.07
TCGA.BH.A0HA.01A.11R.A12P.07            TCGA.BH.A0HA.01A.11R.A12P.07
TCGA.A2.A0YI.01A.31R.A10J.07            TCGA.A2.A0YI.01A.31R.A10J.07
TCGA.A2.A0EX.01A.21R.A034.07            TCGA.A2.A0EX.01A.21R.A034.07

$ tail CV_3.query-results.report.txt.PRANK 
TCGA.BH.A0B9.01A.11R.A056.07    7.20    TCGA.BH.A0B9.01A.11R.A056.07
TCGA.AR.A1AQ.01A.11R.A12P.07    7.19    TCGA.AR.A1AQ.01A.11R.A12P.07
TCGA.B6.A0RT.01A.21R.A084.07    7.17    TCGA.B6.A0RT.01A.21R.A084.07
TCGA.A2.A0D0.01A.11R.A00Z.07    7.17    TCGA.A2.A0D0.01A.11R.A00Z.07
TCGA.A7.A13D.01A.13R.A12P.07    7.16    TCGA.A7.A13D.01A.13R.A12P.07
TCGA.C8.A12K.01A.21R.A115.07    6.93    TCGA.C8.A12K.01A.21R.A115.07
TCGA.A2.A0YM.01A.11R.A109.07    6.87    TCGA.A2.A0YM.01A.11R.A109.07
TCGA.AR.A0U4.01A.11R.A109.07    6.86    TCGA.AR.A0U4.01A.11R.A109.07
TCGA.A8.A07O.01A.11R.A00Z.07    6.84    TCGA.A8.A07O.01A.11R.A00Z.07
TCGA.AN.A0XU.01A.11R.A109.07    6.65    TCGA.AN.A0XU.01A.11R.A109.07
```

<a name="profile"></a>
#### .profile
Feature-level patient data that gets converted into a single similarity network. For example, the dataset contains gene expression and networks are at the level of a pathway, then each pathway (network) would have an input profile of gene expression limited to just the genes for that pathway. These are temporarily created by netDx and converts to similarity networks by calling the GeneMANIA function `ProfileToNetworkDriver`.
Rows are patients and columns are the units grouped together for a network (e.g. genes in a pathway); values are the corresponding patient data.
Tab-delimited file without a header. First column contains ID name for patient. 
In this example, the profile contains information for 3 units (columns) and 3 patients (rows).
```
PatientA   0.25   0.53  -1.224
PatientB    -1.540042       -1.092875       -1.249   
PatientK	-0.290125       -0.212625       -0.825
```

<a name="query"></a>
#### .query
GeneMANIA query file as required by QueryRunner.
File format (as described in the [QueryRunner man page](http://pages.genemania.org/tools/#query-runner)):

Note: the first line is hard-coded by netDx to be 'predictor'. Do not change this unless you understand what you are doing.
```
predictor
query-patient-1 [ \t query-patient-2 ... ]
networks 
related-patient-limit
[ combining-method ]
```
Example:
```
predictor
TCGA.AO.A12G.01A.11R.A10J.07    TCGA.E2.A156.01A.11R.A12D.07    TCGA.C8.A132.01A.31R.A115.07    TCGA.A8.A083.01A.21R.A00Z.07    TCGA.AO.A12E.01A.11R.A10J.07    TCGA.BH.A18S.01A.11R.A12D.07    TCGA.E2.A14Q.01A.11R.A12D.07    TCGA.A2.A0EX.01A.21R.A034.07    TCGA.AO.A12A.01A.21R.A115.07
all
232
automatic
```

<a name="cvscore"></a>
#### pathway_CV_score.txt or pathway_cumTally.txt
This file contains feature-level scores following feature selection. *Note:* The pathway-centric terminology is outdated and will be replaced in future versions of netDx. The results will appear the same no matter what the features are, and they are not limited to pathways. 
Tab-delimited file with a header row. 'PATHWAY_NAME' and 'SCORE'.

```
PATHWAY_NAME    SCORE
CELLULAR_RESPONSES_TO_STRESS_cont.txt   19
NABA_ECM_AFFILIATED_cont.txt    16
COPII__COAT_PROTEIN_2__MEDIATED_VESICLE_TRANSPORT_cont.txt      10
SIGNALING_BY_WNT_cont.txt       13
SIGNALING_BY_TGF-BETA_RECEPTOR_COMPLEX_cont.txt 10
TRANSCRIPTIONAL_ACTIVITY_OF_SMAD2_SMAD3:SMAD4_HETEROTRIMER_cont.txt     10
L1CAM_INTERACTIONS_cont.txt     21
GAP_JUNCTION_ASSEMBLY_cont.txt  10
GAP_JUNCTION_TRAFFICKING_AND_REGULATION_cont.txt        20
```
