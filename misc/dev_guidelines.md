# Developing for netDx

This page has guidelines for developing functions for the netDx software package. Use this as a checklist when writing new functionality. A lot of these follow current best practices in the R software dev community.

## Checklist
1. Document your functions. Make sure you have all the <a href="#tag">basic tags</a>.
2. Make sure the input/output data formats are completely specified.
3. Lines should wrap at 80 char. Set your code editor to highlight column 80/wrap.
4. Make sure there are no outdated/unnecesary input args.
5. Make sure there is no junk code in the function.
6. Your @example may need new data to be added to the package. If so, add that.
7. Rebuild and <a name="testit">test the new functionality</a>. Broken packages add frustration. Build the documentation and the package (`devtools::document()`) and the package (`R CMD INSTALL`). Re-attach the new package (`detach(package:netDx,unload=T);require(netDx)`. 
8. Repeat steps 1-7 till ready. 
9. When ready, generate a pull request.

<a name="tag"></a>
## Basic tags
This package uses the roxygen2 documentation framework to generate Rd files from inline documentation. 
This tutorial is helpful. http://kbroman.org/pkg_primer/pages/docs.html

Required tags for functions:
 - title
 - @details if necessary
 - @params for all. For format of this, see one of the netDx functions.
 - @return output variable
 - @export
 - @examples

Here is example documentation:
```
#' Create patient networks from full matrix of named measurements
#'
#' @details Creates patient similarity networks when full matrices of 
#' data are provided (e.g. gene expression, questionnaire results). To
#' generate networks from sparse data such as CNVs or indels, use 
#' \code{makePSN_RangeSets} instead.
#' The rows of the data matrix (xpr) must be named (nm); one network is 
#' create for each named set (namedSets). There are two options for the 
#' way in which networks are created, depending on the value of
#' \code{writeProfiles}. 
#' 1. writeProfiles=TRUE: GeneMANIA is used to generate interaction networks
#' and sparsify networks. This only works if the desired measure of
#' similarity is network-level Pearson correlation; an example is networks
#' at the level of pathways. In this case, the user does not explicitly 
#' specify a similarity measure and \code{simMetric} is ignored.
#' 2. writeProfiles=FALSE: GeneMANIA is not used to generate interaction
#' networks. Rather, netDx uses \code{simMetric} to create interaction
#' networks. Networks can be sparsified by excluding weak connections 
#' (cutoff). 
#' @param xpr (matrix) rows are measurements, columns are samples. Columns
#' must be named (patient ID)
#' @param nm (character) names for measurements corresponding to row order
#' of \code{xpr}. Must match the names in the named sets specified in
#' \code{nameSets}
#' @param namedSets (list) sets of names to be grouped together. keys are
#' set names, and networks will be named as these. values are character
#' vectors corresponding to groups of names (matching those in \code{nm})
#' that are input to network generation
#' @param outDir (char) path to directory where networks are written
#' @param simMetric (char) measure of similarity. See \code{getSimilarity()}
#' for details
#' @param cutoff (numeric) patients with similarity smaller than this value
#' are not included in the corresponding interaction network
#' @param verbose (logical) print detailed messages
#' @param numCores (integer) number of cores for parallel network generation
#' @param writeProfiles (logical) use GeneMANIA's ProfileToNetworkDriver to
#' create interaction networks. If TRUE, this function writes subsets 
#' of the original data corresponding to networks to file (profiles). 
#' If FALSE, uses  getSimilarity() and writes interaction networks.
#' @param sparsify (logical) sparsify networks by calling sparsifyNets()
#' with default parameters. Only used when writeProfiles=FALSE
#' @param append (logical) if TRUE does not overwrite netDir.
#' @param ... passed to \code{getSimilarity()}
#' @return (char) Basename of files to which networks are written.  
#' Side effect of writing interaction networks in \code{outDir}
#' @examples data(TCGA_mini,pathwayList); 
#' # you may get a warning message that the output directory already
#' # exists; ignore it
#' out <- makePSN_NamedMatrix(xpr,rownames(xpr),pathwayList, 
#' 	".",writeProfiles=TRUE)

#' @export
```

<a name="testit"></a>
## Test package
Have three windows open in your workspace:
1. R session: One to rebuild the doc
2. Terminal: One to reinstall the package. 
3. R session: One to test the reinstalled package.

Say the netDx repo is in: `/Users/Toto/software/netDx`

In 1: 
```
> getwd()
"/Users/Toto/software/netDx/netDx"
> devtools::document()
```
In 2: 
```
$ pwd 
/Users/Toto/software/netDx
$ R CMD INSTALL netDx
```
In 3: 
```
> getwd()
"/Users/Toto/Analysis"
> detach(package:netDx,unload=TRUE);require(netDx) # detach/reattach pkg
> runExample(test_data)
```
Great tutorial on quickly building an R package: https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/
