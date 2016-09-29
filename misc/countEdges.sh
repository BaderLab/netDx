#!/bin/bash

## 10-fold CV, gene expression+CNV
netDir=/home/spai/tmp/TCGA_BRCA/tmp
## 3-way resampling, gene expression only
netDir=/home/spai/tmp/TCGA_BRCA_geneXpr_resample_test/tmp

# count lines in CNV files (_cont.txt)
# get the first field which has num rows
# add num rows (except last line which is total)
wc -l ${netDir}/*_cont.txt | tr -s ' ' | cut -f 2 -d' ' > tmp.txt
sed -i '$ d' tmp.txt
echo "CNV nets\n"
echo "------------------"
echo ""
cat tmp.txt | awk '
	BEGIN {sum=0; maxie=0;minnie=1000000} 
		{
			sum+=$1; 
			if ($1 < minnie) {minnie=$1}
			if ($1 > maxie) {maxie=$1}
		} 
		END {print NR" files; "minnie"-"maxie" edges (mean="sum/NR")"}'
rm tmp.txt

echo "All nets"
echo "------------------"
echo ""
wc -l ${netDir}/INTERACTIONS/*.txt | tr -s ' ' | cut -f 2 -d' ' > tmp.txt
sed -i '$ d' tmp.txt
cat tmp.txt | awk '
	BEGIN {sum=0; maxie=0;minnie=1000000} 
		{
			sum+=$1; 
			if ($1 < minnie) {minnie=$1}
			if ($1 > maxie) {maxie=$1}
		} 
		END {print NR" files; "minnie"-"maxie" edges (mean="sum/NR")"}'
rm tmp.txt
