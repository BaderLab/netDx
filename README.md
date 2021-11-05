# netDx: Network-based patient classifier
[![Docker Build](https://github.com/RealPaiLab/netDx/actions/workflows/push-docker.yml/badge.svg)](https://github.com/RealPaiLab/netDx/actions/workflows/push-docker.yml)
[![R CMD check bioc](https://github.com/RealPaiLab/netDx/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/RealPaiLab/netDx/actions/workflows/check-bioc.yml)

netDx is for biomedical researchers who want to integrate multi-modal patient data to predict outcome or patient subtype. 
netDx builds interpretable machine-learning patient classifiers. Unlike standard machine-learning tools, netDx allows modeling of user-defined biological groups as input features; examples include pathways and co-regulated elements. In addition to patient classification, top-scoring features provide mechanistic insight, helping drive hypothesis generation for downstream experiments. netDx currently provides native support for pathway-level features but can be generalized to any user-defined data type and grouping.

## Install

### Get Dockerized version
Get the latest stable version of netDx with all dependencies installed in a Docker container. 

Connect your code and data to the container by mapping volumes from your machine to Docker and run RStudio in a web browser. See example below.

These steps were tested on a 2020 Mac Book Pro with 32 Gb RAM.

1) Install Docker, and under Resources for Docker Desktop's settings set Memory to 5 GB.
2) In terminal, pull the netDx image from Docker:
	`docker pull realpailab/netdx`
3) To run the image as a container,  save the following code as a bash script named `startDocker.sh`:
```
#!/bin/bash
	
docker run -d --rm \
    -p 8787:8787 \
    -e PASSWORD=netdx \
    -v /path/to/software/dir:/software \
    -v /path/to/data/dir:/data \
    --name [containerName] \
    realpailab/netdx
```
4) Modify the bash script so that it's executable and run it:
```
chmod u+x startDocker.sh
./startDocker.sh
```
5) Access RStudio by going to your web browser and typing in `localhost:8787`, and sign in with the username `rstudio` and password `netdx`.
6) You can then run the vignettes in R:
```
setwd("/home/rstudio/vignettes")
rmarkdown::render("RunPredictorWorkflow.Rmd")
```


### Install in Rstudio


1. If you don't already have it, install [base R](https://www.r-project.org/) and [RStudio Desktop](https://www.rstudio.com/products/rstudio/download/

2. Install Java runtime for your OS.
**Note:** netDx requires *Java 16 or earlier*; it will not work with Java 17 or higher.
If you have Java 16 or earlier, go to step 3.

**Unix:** Find and install the latest JDK and JRE libraries by running this on terminal:
```
apt-cache search jdk
sudo apt-get install openjdk-11-jre openjdk-11-jdk
```
This operation requires system administrator privileges. You are done, move to Step 3.

**OS X:** 
[Download Java 16.0.1 (build 16.0.1+9) for 64-bit Mac](https://download.java.net/java/GA/jdk16.0.1/7147401fd7354114ac51ef3e1328291f/9/GPL/openjdk-16.0.1_osx-x64_bin.tar.gz) and run the installer.

Install it by opening terminal and typing the following:
```
sudo mv openjdk-16.0.1_osx-x64_bin.tar.gz /Library/Java/JavaVirtualMachines/
cd /Library/Java/JavaVirtualMachines/
sudo tar -xzf openjdk-16.0.1_osx-x64_bin.tar.gz
sudo rm openjdk-16.0.1_osx-x64_bin.tar.gz
```
The above operation requires `sudo` privileges (system admin privileges) on the machine.

Now in Rstudio, change the value of `JAVA_HOME` to point it to this version of Java:
```
Sys.setenv(JAVA_HOME="/Library/Java/JavaVirtualMachines/jdk-16.0.1.jdk/Contents/Home")
```

Check to make sure R is using Java 16.0.1:
```
system2("java","--version")
```

If you see this, you're set!
```
openjdk 16.0.1 2021-04-20
OpenJDK Runtime Environment (build 16.0.1+9-24)
OpenJDK 64-Bit Server VM (build 16.0.1+9-24, mixed mode, sharing)
```

3. In R, [install BioConductor](https://www.bioconductor.org/install/) if needed. 
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
```

4. Install netDx.
```
BiocManager::install("netDx",dependencies=TRUE)
```


---

### Main repo for netDx dev work as of Sep 2021.

netDx is a general-purpose algorithm for building patient classifiers by using patient similarity networks as features. It excels at interpretability and handling missing data. It also allows custom grouping rules for features, notably grouping genes into pathways. It integrates with RCy3 for network visualization of predictive pathways.

As of February 2020, netDx is available via the BioConductor repository. 
Visit http://bioconductor.org/packages/release/bioc/html/netDx.html to install the package and see worked examples.

Contact Shraddha Pai at shraddha.pai@utoronto.ca in case of questions.


References: 

1. Pai S, Hui S, Isserlin R, Shah MA, Kaka H and GD Bader (2019). netDx: Interpretable patient classification using patient similarity networks. *Mol Sys Biol*. 15: e8497. [Read the paper here](https://www.embopress.org/doi/full/10.15252/msb.20188497).
2. Pai S, Weber P, Isserlin R, Kaka H, Hui S, Shah MA, Giudice L, Giugno R, Nøhr AK, Baumbach J, GD Bader (2021). netDx: Software for building interpretable patient classifiers by multi-'omic data integration using patient similarity networks. *F1000 Research*. 9:1239.
