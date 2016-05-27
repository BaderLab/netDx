# netDx
Patient classifier based on patient similarity networks

This repo contains two R packages, each of which need to be separately installed:

1. `netDx/`: Software implementing the netDx method
2. `netDx.examples/`: Data for use cases

The `examples/` folder contains R code that should just run once both `netDx/` and `netDx.examples/` are installed.

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
From command-line, download the git repo for these packages and install them:

```
$ git clone git@github.com:BaderLab/netDx.git
$ cd netDx/
$ R CMD INSTALL netDx
$ R CMD INSTALL netDx.examples
```

## Run examples
Each example runs as a standalone `.R` function. To run these, you need to have both `netDx` and `netDx.examples` installed. From the directory where the `netDx/` repo was downloaded:
```
$ cd netDx/examples/
$ R
R version 3.2.4 (2016-03-10) -- "Very Secure Dishes"
Copyright (C) 2016 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin13.4.0 (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

  Natural language support but running in an English locale

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> source("TCGA_BRCA.R")
```
