# SEAPODYM R tools: SEAMOVE library

R library to extract SEAPODYM movement fluxes (biomass flow rate) from donor to recepient region as a result of seapodym\_fluxes simulation with specific regional structure, and to convert them to the seasonal movement probabilities to be used in Multifan-CL. The regional and age structure of Multifan-CL should be defined in the SEAPODYM configuration (XML parfile) for the simulation.

## About

Package: seamove
Version: 0.0
Title: Reading and visualizing SEAPODYM regional fluxes as well as exporting them in the Multifan-CL format.
Author: Inna Senina <innas@spc.int>
Maintainer: Inna Senina <innas@spc.int>
Depends: R (>= 1.0.0)
Description: This small package provides functions to read outputs of seapodym\_fluxes application, from a single or multiple simulations, plot them and to export them in the format used by Multifan-CL, i.e. as seasonal movement probabilities for a given age and regional structure.

## Installation

Method #1. 

1. In command line outside of downloaded package directory (Linux only)

   CMD build seamove

2. From pkg.tar.gz file, in command line (prompt in Windows)

   R CMD INSTALL seamove\_0.0.tar.gz

   Note, in Windows install the Rtools package to enable the installation from tar.gz

Method #2. With R devtools package, in R environment (no need to build first)

   library(devtools)
   devtools::install\_github("[path]/seapodym/rtools/seamove")

## Algorithm


## Examples 

Simulations of seapodym\_fluxes were run with the following regional structure specified in the XML parfile.

![Albacore assessment regions](/images/alb-regions.png)

The simulation directory `ALB.DIR` contains the folder `output` with alb\_FluxesRegion\_age[a].txt.

  `dir.in <- paste0(ALB.DIR,"/output/")`
  `dir.out<-paste0(ALB.DIR,"/output-mfcl/")`

Let's run the main function handling all necessary steps to extract, convert the units and plot the resulting movement probabilities as computed for the Multifan-CL model:

  `sea2mfcl.movement(dir.in,dir.out,1979:2010,age.in=5:148,age.plus=FALSE,aggregate.age=3,nbr.in=4)`

In this example, the monthly age classes starting from age class 5 (sixth age class with mean age 5.5 months, the information provided in the outputs) till age class 148 of SEAPODYM (the value of age.plus=FALSE indicates that this is not an A+ class) will be aggregated to quarterly resolution (parameter `aggregate.age=3`), resulting in 48 age classes. The routine will write the movement-matrices-for-MFCL.txt file as well as plot the results as shown on the figure below for the movement from region 3 to region 2.

![](/images/movement_probability_r3-to-r2.png)



