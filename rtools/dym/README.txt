SEAPODYM R tools: DYM library

-----------------------------------------------------------
Package: dym
Version: 0.0
Title: Read and write DYM files
Author: Inna Senina <innas@spc.int>
Maintainer: Inna Senina <innas@spc.int>
Depends: R (>= 1.0.0)
Description: This small package provides functions to read, 
             write and extract the subsets of data from DYM 
	     files that are the inbuilt binary files of the 
	     SEAPODYM model.
-----------------------------------------------------------

-----------------------------------------------------------
Build:

1. In command line (linux, for Windows only see method #2)

   CMD build dym  

Installation:   

1. From package.tar.gz file, in command line (prompt in Windows)

   R CMD INSTALL dym_0.0.tar.gz

   Note, in Windows install the Rtools package to enable the installation from tar.gz file

2. With R devtools package, in R environment

   library(devtools)
   devtools::install_github("PacificCommunity/seapodym-codebase/tree/master/rtools/dym")
  



