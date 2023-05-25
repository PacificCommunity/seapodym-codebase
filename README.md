# SEAPODYM project

Source code of the SEAPODYM numerical modelling framework and analyses tools, OFP-FEMA. C++, R, Linux.

*Version 4.1*

*Pre-released January 21, 2023*

## About

SEAPODYM (Spatial Ecosystem And POpulation DYnamics Model) has been initiated in the mid 1990s by the Oceanic Fisheries Programme of the Pacific Community (SPC), and its developments were continuously conducted since then at the University of Hawaii (2004-2007), at Collecte Localisation Satellites in Toulouse (2006-2020) and at SPC (2020-now). The objective of this project is to propose a new scientific tool for biomass estimation at spatial scales smaller than operating stock assessment models, while having quantitative skill of predicting fish stocks spatiotemporal dynamics under the influence of both fishing and environmental variability. 

The codebase implements the SEAPODYM modelling framework, which includes several models and tools, allowing running the full fish population dynamics model with pre-configured parameters, estimating this model parameters by integrating geo-referenced fisheries and tagging data, computing and estimating parameters of spawning and feeding habitats, simulating movement dynamics of tagged fish, and evaluating biomass flow rates between designed regions. The model can be parametrised through integration of the commercial fisheries data and the scientific tagging data, and solving optimization problem formulated using the maximum likelihood estimation approach. To ensure the optimization algorithm efficiency, implementation of the model adjoint code allows an exact, analytical evaluation of the likelihood gradient. Additional tools integrated within C++ code include local and global sensitivity analysis, and Hessian approximation. The offline R tools, developed for various model analyses and the IO files manipulations, are included in the repository to facilitate the use of the model.

## Prerequisites 

SEAPODYM runs on a 64-bit computer and on Linux operating system only. As all highly dimensional numerical applications, exploitation of SEAPODYM puts requirements on available physical memory and CPU, which depend on model spatiotemporal resolutions. For example, ready-to-use configurations of SEAPODYM models at 1 degree square and 1 month resolution can run on a laptop with 8GB of RAM and a CPU with 2GHz frequency. For parameter estimations and for higher resolution simulations, it is advised to use higher performance computers with at least 32 Gb of RAM. Note, current version of SEAPODYM uses only one CPU per instance. 

## Installation

1. **Precompiled binaries**

  Supported operating systems are 

  * Ubuntu 
  * Red Hat Linux

  Note, the compiler GNU C++ should be installed. 

  Download the precompiled binary release corresponding to your Linux system and gcc versions from [Releases](https://github.com/PacificCommunity/seapodym-codebase/releases/tag/seapodym-4.1/) directory. Extract contents of the zip file to a folder _~/seapodym/_.

  In the **Terminal** window, go inside the extracted *bin* folder and type

  [~/seapodym/bin/]$ seapodym -h

  This command should display the usage instruction and running options of SEAPODYM.

2. **Building from source**

  Source code compilation requires the GNU C++ compiler installed on the computer. 

  In addition, the following two libraries should be installed:

  * libxml2 
  * AUTODIF 

  The libxml2 library is used to read and write all application parameters in the XML file. If not installed already, you can install it with the command:

    [~]$ sudo apt-get install libxml2-dev

  choosing the package libxml2-dev, which contains header files that are required to compile SEAPODYM applications.   

  The AUTODIF libraries provide an array language extension to C ++ and can be installed as a part of the ADModel Builder software. Visit the [ADMB project](https://github.com/admb-project/admb) page and follow the instructions for [quick installation](https://github.com/admb-project/admb/blob/main/docs/install/QuickStartUnix.md) or [building from source](https://github.com/admb-project/admb/blob/main/docs/install/BuildingSourceUnix.md) of the latest release of the ADMB software. 

  Once ADMB software is installed, in the .bashrc file declare environment variable pointing to it, as follows:

    
    export ADMB_HOME=/your-path-to-admb-folder/
    

  In the case when shared (dynamic) libraries are to be used, add the following line in the .bashrc:


    export LD_LIBRARY_PATH=$ADMB_HOME/lib:$LD_LIBRARY_PATH


  Finally, create the following soft link for the optimized AUTODIF library in the _lib_ folder that was compiled with the gcc version installed on your machine. If the directory containing ADMB package is _~/admb/_, the static library files are located in _~/admb/lib/_. The one that is required by _seapodym_ is called \*_libadmbo_\*.a, for example _libadmbo-x86_64-linux-g++11.a_. To create a symbolic link do


    [~/] cd ~/admb/lib

    [~/admb/lib/]$ln -s libadmbo-x86_64-linux-g++11.a libadmbo.a


  Once the required libraries have been installed and configured, proceed to installing SEAPODYM. The source code can be downloaded either from the [latest release](https://github.com/PacificCommunity/seapodym-codebase/releases/seapodym-4.1/), or, alternatively, checking out the github repository:

    
    [~]$ git clone https://github.com/PacificCommunity/seapodym-codebase.git
    
  
  If downloaded a release, unpack the archive. Go into the folder with the *src* directory and Makefiles, which is presumably named _~/seapodym/_. Compile one of the SEAPODYM applications by typing the command
    
    [~/seapodym/]$ make
    
  If the compilation is successful, it will create the binary file _seapodym_ in directory _bin_. This application runs the model in a simulation mode only. Executing
   
    [~/seapodym/bin/]$ seapodym -h
  
  should display the usage instruction and running options. To compile the sub-model of species habitats type

    [~/seapodym/]$ make -f Makefile.hab

  The binary file *seapodym\_habitats* should be created. Type

    [~/seapodym/bin/]$ seapodym_habitats -h

  to see the running options of this SEAPODYM sub-model. To compile all model applications at once, run the following batch file inside folder containing Makefiles

    [~/seapodym/bin/]$ . all.bat

  For convenience, the path to the _~/seapodym/bin/_ directory can be added to variable PATH in the user's _~/.bashrc_. 


## Running SEAPODYM models
    
  Example #1. **Habitat**

  Go to the *example-configs* directory and test that the binaries execute nominally and the models generate expected outputs on your computer by executing a small habitat model example
    
    [~]$ cd ~/seapodym/examples-configs/habitat    

  **Simulation** run

  To execute habitat model in simulation mode type the command with option *-s*
    
    [~/seapodym/examples-configs/habitat]$ seapodym\_habitats -s habitat.xml

  The simulation log can be compared with the one provided in the sim.out file. 
  
  **Optimization** run

  If option *-s* is omited, the application _seapodym\_habitats_ will start optimization. Note, due to much higher CPU and RAM demands compared to simulation mode, running optimization experiments is advised only after familializing yourself with the model and method basics, described in [Model Reference Manual](docs/manual/Seapodym_user_manual.pdf). 

  The optimization in _seapodym\_habitats_ aims to estimate habitat parameters by fitting the modelled habitat to an a priori known (e.g., previously estimated) habitat field. This can be convenient when only habitat parameters need to be (re-)estimated, e.g., in the twin experiments, or once the ocean forcing fields or other external parameters had been modified. For this example, "observations" were generated by running _seapodym\_habitats_ simulation with parameters in file _habitat\_input.xml_. The 3D (time and two-dimensional space) habitat field is written in binary file _msp\_spawning\_habitat\_input.dym_ (currently stored in directory _output_), which can be read with help of [dym](https://github.com/PacificCommunity/seapodym-codebase/tree/master/rtools/dym) R library. To test this application in optimization mode, run it starting with non-optimal parameters

    [~/seapodym/examples-configs/habitat]$ seapodym\_habitats habitat.xml

  The output should be the same as in the log file _opt.out_. The new parameters written in _newparfile.xml_ once the optimization converged to a minimum should be very close to that in habitat\_input.xml.  


  Example #2. **Skipjack**

  This is a more comprehensive full model example, the pre-configured model of skipjack tuna, based on the [parameter estimation approach with fisheries and tagging data](https://cdnsciencepub.com/doi/full/10.1139/cjfas-2018-0470). To run this model, the corresponding forcing directory needs to be downloaded from the public [data repository](https://osf.io/h8u93) on the OSF platform.

  Unzip and place the forcing files into a local directory without modifying the folder structure.

  Now open the skipjack\_F0.xml parfile in the preferred text editor and modify the ${SEAPODYM\_HOME}, so to provide the full absolute path to the unzipped *data* folder.  
    
    [~]$ cd example-configs/skipjack/

  **Simulation** run    
    
  To run this optimal model in simulation mode type the following command

    [~/seapodym/example-configs/skipjack/]$ seapodym -s skipjack_F0.xml
    
  The simulation log can be compared with the one provided in sim\_F0.out file. The binary outputs will be written in folder _output/output\_F0\_/_.

  Note that running this simulation with fishing requires fisheries data, which are not public. To inquire for the access to the data, please email <ofpdatarequest@spc.int>. Your request will be examined based on the eligibility of your organisation and your project. Then, the simulation with fishing can be run with the XML configuration file skipjack.xml, and the screen log compared with provided sim.out file.

  **Optimization** run

  It is possible to run optimization without fisheries data using _seapodym\_densities_ application. Make sure it has been built. If not run the batch file to compile all five applications (see *Building from source* instructions above). Similarly to _seapodym\_habitats_, this application uses model outputs as "observations". Using this example configuration, make a copy of _skipjack\_F0.xml_ parfile

    [~/seapodym/example-configs/skipjack/]$ cp skipjack_F0.xml initparfile.xml

  Now reduce the simulation time period to five years by modifying the _year_ attribute in _save_last_date_ node from 2010 to 1983. First run the simulation with this parfile

    [~/seapodym/example-configs/skipjack/]$ seapodym -s initparfile.xml

  Make the input file with pseudo-observations by copying the binary output file _skj\_totbm.dym_ 
    
    [~/seapodym/example-configs/skipjack/]$ cp output/output_F0/skj_totbm.dym output/output_F0/skj_density_input.dym

  Now you can run optimization

  [~/seapodym/example-configs/skipjack/]$ seapodym\_densities initparfile.xml

  and compare its outputs with that provided in log file _opt.out_.  
  

  Example #3. **Albacore**

  Similarly to skipjack example, this is the pre-configured example to run the albacore tuna model, [parameterized and validated with fisheries data](https://www.sciencedirect.com/science/article/pii/S0967064519301511). The configuration file will allow running the population dynamics model of the south Pacific stock. To run this model, download the forcing directory from the [OSF data repository](https://osf.io/j53hc).

  Unzip and place the forcing files into a local directory without modifying the folder structure.

  Open the albacore\_F0.xml parfile and modify the ${SEAPODYM\_HOME}, so to provide the full absolute path to the unzipped *data* folder. 
    
    [~]$ cd example-configs/albacore/
    
  and run 

    [~/seapodym/example-configs/albacore/]$ seapodym -s albacore_F0.xml
    
  The simulation log can be compared with the one provided in sim\_F0.out file. Note that running this simulation with fishing requires fisheries data, which are not public. To inquire for the access to the data, please email <ofpdatarequest@spc.int>. Your request will be examined based on the eligibility of your organisation and your project. Then, the simulation with fishing can be run with the XML configuration file albacore.xml, and the screen log compared with provided sim.out file.

     
## Documentation
For further information on how to use SEAPODYM for simulation studies and model developments see the [Model Reference Manual](docs/manual/Seapodym_user_manual.pdf). Doxygen generated [Code Documentation](docs/code-dox/codedoc_seapodym.pdf) can be consulted for a quick overview of the C++ numerical model code. 

## License
SEAPODYM is an open source project. The model code, associated tools and documentations are provided under the general terms of the three-clause BSD [LICENSE](LICENSE.md).
