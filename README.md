# netDx
netDx is a method for building patient classifiers based on patient similarity networks.

This repo contains two R packages, each of which need to be separately installed:

1. `netDx/`: Software implementing the netDx method
2. `netDx.examples/`: Data for use cases

The `examples/` folder contains R code that should just run once both `netDx/` and `netDx.examples/` are installed.

For more on this project, visit **[http://netdx.org]**

## Installation

### Prerequisites
You must have Java, Python and R installed. Within R, you must have BioConductor installed. This section helps you figure out which
of these you need to install. If you already have all these, skip to the next section.

#### Java (1.8+ recommended, but will probably work on 1.6+)

At command line, run `java --version`. You should see output like this:
```
java version "1.8.0_31"
Java(TM) SE Runtime Environment (build 1.8.0_31-b13)
Java HotSpot(TM) 64-Bit Server VM (build 25.31-b07, mixed mode)
```
If you don't see this kind of output, you may need to first [install java](https://java.com/en/).

#### Python (2.7 recommended)

At command line, run `python --version`. You should see output like this:

```
Python 2.7.9 :: Anaconda 2.1.0 (x86_64)
```

#### R (3.2+ recommended)
At command line, run `R --version`. You should see output like this:
```
R version 3.2.4 (2016-03-10) -- "Very Secure Dishes"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin13.4.0 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under the terms of the
GNU General Public License versions 2 or 3.
For more information about these matters see
http://www.gnu.org/licenses/.
```

If not, [install R](https://www.r-project.org/).

#### BioConductor (2.3+ recommended)
[BioConductor](http://bioconductor.org/) is a system of R objects and software specifically for biological applications.
Once you have R installed, install the `biobase` and `GenomicRanges` packages via the BioConductor installer:
```
$ R
> source("http://bioconductor.org/biocLite.R")
> biocLite(c("Biobase","GenomicRanges"))
```
Say `yes` to all dependencies that need to be installed.

## Install `netDx` and `netDx.examples`
This section assumes you have Java, Python, R and Bioconductor installed. From command-line, download the git repo for these packages and install them. In the code below, output from intermediate steps is omitted for clarity.

*Note: For now, R package dependencies must be separately installed using the install.packages() call as shown below. netDx will be submitted to CRAN following publication; thereafter, dependencies can be automatically installed with the call to install netDx.*

```
$ git clone git@github.com:BaderLab/netDx.git
$ cd netDx/
$ R
> install.packages(c("bigmemory","foreach","combinat","doParallel","ROCR","pracma","RColorBrewer","reshape2"))
> install.packages("netDx",type="source",repos=NULL)
> install.packages("netDx.examples",type="source",repos=NULL)
```

## Run examples
Each vignette is in Sweave format (`.Rnw`) . To run these, you need to have both `netDx` and `netDx.examples` installed, and `knitr` to reproduce the results. From the directory where the `netDx/` repo was downloaded, here is code that runs the breast cancer example using data from The Cancer Genome Atlas (Ref 1):
```
$ cd netDx/examples/
$ R
> require(knitr)
> knit2pdf("Medulloblastoma.Rnw")
```
This should generate `Medulloblastoma.pdf` in the `examples/` directory. 
Alternatively, if you have [Rstudio](https://www.rstudio.com/home/) installed (highly recommended), you should be able to open the `Rnw` file and click `Compile PDF`. If you are getting an error saying that R cannot find `/usr/bin/texti2dvi`, install the `texinfo` package.

### References
1. The Cancer Genome Atlas Network. (2012). Comprehensive molecular portraits of human breast tumours. Nature. 490:61.
