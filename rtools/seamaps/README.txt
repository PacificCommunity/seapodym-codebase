SEAPODYM R tools: SEAPODYM-MAPS library

------------------------------------------------------------
Package: seamaps
Version: 0.1
Title: Nice maps with land, coastline and grid, which can be 
       combined with any geo-referenced data (oceanographic,  
       biological, catch, etc.).
Author: Inna Senina <innas@spc.int>
Maintainer: Inna Senina <innas@spc.int>
Depends: R (>= 1.0.0)
Description: This small package provides functions to plot 
             nice maps, color palettes and utilities to 
	     prepare SEAPODYM data for plotting.
-----------------------------------------------------------

-----------------------------------------------------------
Build (Linux only, for Windows see method #2):

1. In command line outside of downloaded package directory

   CMD build seamaps

Installation:   

1. From pkg.tar.gz file, in command line (prompt in Windows)

   R CMD INSTALL seamaps_0.0.tar.gz

   Note, in Windows install the Rtools package to enable the 
   installation from tar.gz

2. With R devtools package, in R environment (no need to 
   build first)

   library(devtools)
   devtools::install_github("[path]/seapodym/rtools/seamaps")
  



