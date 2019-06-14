#!/bin/bash

# Script to automate dependency install for netDx

echo "* Checking if Java installed ..."
if java -version 2>&1 > /dev/null |  grep -q "java version" ; then
  echo -e "\tdone."
else {
  echo -e "*** ERROR: Java not found; install (https://www.java.com/en/download/) or add to path"
  exit 0;
 }
fi

echo "* Checking if Python installed ..."
if [[ $(python --version 2>&1)  ]]
    then
        echo -e "\tdone"
    else {
        echo -e "*** ERROR: Python not found; install (https://www.python.org/downloads/) or add to path"
	exit 0;
	}
fi

echo "* Checking if R installed ..."
if R --version | grep -q "R version" ;  
    then
			 ver=`R --version | grep "R version" | cut -f 3 -d " "`
			echo -e "\tversion found: $ver"
	 	   ver1=`echo $ver | cut -f1 -d"."`
		   ver2=`echo $ver | cut -f2 -d"."`
			if [ $ver1 -ge 3 ] &&  [ $ver2 -ge 6 ]; then
        echo -e "\tdone"
		  else {
				echo ""
				echo -e "\t*** ERROR: Version 3.6+ of R required. Install from https://cran.r-project.org/, or add to path"
			  exit 0
		}
			fi
    else {
				echo -e "\t*** ERROR: R not found. Install R 3.6+ from https://cran.r-project.org/, or add to path"
	exit 0;
	}
fi

# install R packages
echo "* Installing R dependencies"
echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile

declare -a PKGS=( devtools curl bigmemory foreach combinat doParallel ROCR pracma RColorBrewer reshape2 ggplot2 caroline rmarkdown igraph glmnet );
for p in ${PKGS[@]};do 
	echo -e "\t* Checking for $p"
	Rscript -e "if(!requireNamespace(\"$p\",quietly=TRUE)){ install.packages(\"$p\")}"
done

echo "* Installing BioConductor if required"
Rscript -e 'if (!requireNamespace("BiocManager", quietly = TRUE)){install.packages("BiocManager")}'

echo "* Installing BioConductor dependencies if required"
declare -a PKGS=( GenomicRanges RCy3 );
for p in ${PKGS[@]};do 
	echo -e "\t* Checking for $p"
	Rscript -e "if(!requireNamespace(\"$p\",quietly=TRUE)){ BiocManager::install(\"$p\")}"
done

echo "* Checking if pandoc installed (needed to run tutorials) ..."
if pandoc -v | grep -q "^pandoc " ;  
    then
			 ver=`pandoc -v | grep "^pandoc " | cut -f 2 -d " "`
			echo -e "\tversion found: $ver"
	 	   ver1=`echo $ver | cut -f1 -d"."`
		   ver2=`echo $ver | cut -f2 -d"."`
			if [ $ver1 -ge 2 ] ; then
        echo -e "\tdone"
		  else {
				echo ""
				echo -e "\t*** Version 1.12.3+ of pandoc not found! Installing..."
				curl -L https://github.com/jgm/pandoc/releases/download/2.7.2/pandoc-2.7.2-macOS.pkg -o pandoc.pkg
				sudo installer -pkg pandoc.pkg -target /
		}
			fi
   else {
				echo -e "\t*** Version 1.12.3+ of pandoc not found! Installing..."
				curl -L https://github.com/jgm/pandoc/releases/download/2.7.2/pandoc-2.7.2-macOS.pkg -o pandoc.pkg
				sudo installer -pkg pandoc.pkg -target /
	}
fi

cd ..
echo "* Installing netDx"
R CMD INSTALL netDx
