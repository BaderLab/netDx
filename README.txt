# netDx
netDx is an algorithm for building patient classifiers by using patient similarity networks as features.

This repo contains two R packages, each of which need to be separately installed:

1. `netDx/`: Software implementing the netDx method
2. `netDx.examples/`: Data for use cases

The `examples/` folder contains R code that should just run once both `netDx/` and `netDx.examples/` are installed.

For more information and FAQ, visit **http://netdx.org**

* [Install netDx](#install-netdx)
* [Test functionality](#test-functionality)
  * [Known issues with compiling pdfs](#known-issues-with-compiling-pdfs)
* [Run breastcancer LumA example](#run-breastcancer-luma-example)


## Install netDx

### Prerequisites
**netDx has been tested on Mac OS/X and on Linux systems. For now we recommend you run netDx on these operating systems.** Future versions of netDx will have Windows support.

You must have Java, Python and R installed. Within R, you must have BioConductor installed. This section helps you figure out which
of these you need to install. If you already have all these, skip to the next section.

#### Java (1.8+ recommended, but will probably work on 1.6+)
netDx uses the GeneMANIA algorithm to integrate patient networks and recommend patients by similarity (Mostafavi and Morris (2008). *Genome Biol* 9:Suppl 1). GeneMANIA is currently implemented in Java, making this interpreter a requirement for netDx. 

At command line, run `java --version`. You should see output like this:
```
java version "1.8.0_31"
Java(TM) SE Runtime Environment (build 1.8.0_31-b13)
Java HotSpot(TM) 64-Bit Server VM (build 25.31-b07, mixed mode)
```
If you don't see this kind of output, you may need to first [install java](https://java.com/en/).

#### Python (2.7 recommended)
netDx uses legacy Python scripts in creating the GeneMANIA database, so for now a Python interpreter is required to run netDx. Future versions of netDx will not have this requirement.

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
$ cd netDx-master/
$ R
> install.packages(c("bigmemory","foreach","combinat","doParallel","ROCR","pracma","RColorBrewer","reshape2"))
> install.packages("netDx",type="source",repos=NULL)
> install.packages("netDx.examples",type="source",repos=NULL)
> install.packages("knitr") # needed to run examples
```

## Test functionality
Run the medulloblastoma vignette to make sure the netDx pipeline works from end to end.
Each vignette is in Sweave format (`.Rnw`) . To run these, you need to have both `netDx` and `netDx.examples` installed. You will also need to install the R package `knitr` to compile the Sweave file.  If you have [Rstudio](https://www.rstudio.com/home/) installed (highly recommended), you should be able to open the `Rnw` file and click `Compile PDF`. Alternately, you may run the vignette through an interactive R session:

```
$ cd netDx/examples/
$ R
> require(knitr)
> knit2pdf("Medulloblastoma.Rnw")
```
This should generate `Medulloblastoma.pdf` in the `examples/` directory. 

## Known issues with compiling pdfs

#### (Linux) 
If you are getting an error saying that R cannot find `/usr/bin/texti2dvi`, install the `texinfo` package.

#### (OS/X): "`pdfLaTex` not found" error
When compiling the pdf, you may get a message saying that `pdfLaTex` is not installed. We have had one such report on OS/X and it is known to occur after an upgrade to OS X Mavericks. The following steps resolved the issue (paraphrased from [this post](http://stackoverflow.com/questions/22081991/rmarkdown-pandoc-pdflatex-not-found)):
Step 1. Check that `/usr/texbin` exists. In terminal, type:
```
cd /usr/texbin
```
Step 2. If you get a "No such file or directory" message, create a symbolic link to your installation's `texbin` file. It may be in `/Library/TeX/Distributions/.DefaultTeX/Contents/Programs/texbin`.
In terminal type:
```
ln -s /Library/TeX/Distributions/.DefaultTeX/Contents/Programs/texbin /usr/texbin
```
Step 3. Now, in the terminal check the value of `echo $PATH`. Make sure that `/usr/texbin` is present. If it isn't present, then you need to add `/usr/texbin` to your PATH variable. This can be done by updating the `PATH` variable in `~/.bashrc`.
However, if you find yourself having to mess with the PATH variable, try reinstalling the [MacTex](http://tug.org/mactex/) package.

## Run BreastCancer LumA example
This vignette is presented in the netDx manuscript. Here we start with 348 primary tumours from the Cancer Genome Atlas, and build a predictor for Luminal A subtype classification (The Cancer Genome Atlas (2012). *Nature.* **490**:61-70).  This example illustrates  feature selection using a simple design in which networks are scored out of 10 based on a single round of 10-fold cross validation. On a MacBook Air laptop (late 2014), this vignette takes ~1.5 hours to run to completion. You may speed it up by running it on a machine with more processors and changing the `numCores` variable in the vignette. We do not recommend running it on a machine with less than 8Gb RAM.

```
$ cd netDx/examples/
$ R
> require(knitr)
> knit2pdf("BreastCancer.Rnw")
```
