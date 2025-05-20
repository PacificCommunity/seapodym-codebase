# SEAPODYM project: branch DEV-SKJ

Branch **dev-skjmodel** is the SEAPODYM source code branch used for the SKJ reference model revision under JRA55 forcing. 

*Pre-version 5.x*

*Evolved, merged and revised April 2024 -January, 2025*

## About

Branch **development** is the evolution of **master**, **dev_larvae_dynpop**, and *dev-skjmodel* branches. In addition to the larval and other early life data integrated into parameter estimation with the MLE approach (**dev_larvae_dynpop** branch) and early life model revision (*dev-skjmodel* branch), it enables the estimation of catchability slopes. Read (4), (5) and (8) below to know how to adapt the parfile for this code version. 

Compared to the **master** and **dev_larvae_dynpop**, **development** branch includes the following **major** (either significantly impacting model solutions and/or involving substantial changes of adjoint functions) changes and improvements:

1) Normalization of *f_prey* function of the spawning habitat, ensuring estimation of the functional relationship in (0,1) range. This change improves estimation of recruitment parameters and isn't optional. Besides, the old models cannot be reparameterized so to provide exactly the same output with the new function. For further details, see functions *Hs_comp_elem* and *Hs_comp* in *spawning_habitat.cpp*.   

2) Enabling arbitrary number of forage groups in optimization mode (previously working only in simulation mode).

3) Stock likelihood term in the form of penalty, allowing setting the maximal stock size (previously it was the best value) which cannot be exceeded.  

4) Model of early larvae (ELM), with time splitting of the first month of fish life and different processes influencing mortality rate during each time-age interval - SST-driven mortality during egg, yolk-sac and early development stage and HS-driven mortality for the rest of the larval stage. Early larvae are then used in the likelihood with larval sampling data. To activate this model, added a set of new parameters in the block *early_larvae_model*. Also, to be able to use it, it is necessary to enter parameters of empirical egg survival function depending on SST. For further details, see function *M_early_larvae* in *mortality_sp.cpp*. The ELM is optional, so in case, if this model is unnecessary, e.g., no larval data available, it can be desactivated. See <early_larvae_model> in **skipjack_1deg.xml** in *examples*, in section **EARLY-LIFE DATA LIKELIHOOD**.

5) Changes in mortality function: i) new parameter *M_larvae_range* added with the same role as *M_mean_range*, but for the second larval stage only, ii) Rage term that allows variable response to environmental conditions expressed as Hj or Ha has parameters in the parfile (parameters *M_max_range_age* and *M_max_range_slope*), and iii) removing arbitrary value Hval = 0.5 for mortality variability, instead suggesting (by default) using the maximal value H = 1 for minimal natural mortality, and all values below 1 gradually increasing these M-at-age values, giving maximal mortality rates in H = 0. It is however possible to set other values as habitat threshold for mortality change. For example, this must be done in case of running the old reference models with this code, then add the line <Mvar_Hval value="0.5"/> to the parfile, section *Fixed parameters*.

6) Changes in recruitment function: now it can be set up to get the number of newborns depending on adult function only. Set new parameter *spawning_in_hs = 0* to swith to this option. This is a necessary part of the new ELM, in which case this parameter is ignored, however it can be useful for all models in general. Currently, the double effect of HS on larval abundance (at spawning and through the variable mortality, which in addition can both increase and decrease mortality through *M_mean_range*), although at different time steps, likely makes the estimation of demographic parameters more difficult.   

7) Added eF scalers to averaging ocean currents according to the time spent in the layer, function *Average_currents*. This is optional and can be set up via flag *scale_forage_ave_currents*. 

8) The catchability equation is revised. Instead of 'dyn' parameters, with the line slope defined as q x dyn, it is now written in a classical form with slope parameter (q\_slope\_fishery in the code) and intercept being the initial catchability at time 0. Note, the previous parfiles still can be used with this branch, but this functionality is preserved for simulation only to enable a smooth transition to the new definition. For all runs, requiring optimizations, it is necessary to use the new function parameterization as follows:

	<P1 skj="0.001"> <!-- catchability at t=0 -->
			<variable min="0.0" max="1" use="true"/> <!-- boundaries for catchability at t=0 -->
			<slope skj="0.0" min="0.0" max="1" use="true"/> <!-- line slope and its boundaries -->
	</P1>

Other **minor** changes (same as **dev_larvae_dynpop** branch) include: 

1) In the larval data predictions, additional multiplier term linking plankton net catchability to mixed-layer depth can be used. Flag *q_mld_larvae* activates its use. For the moment, the function is hard-coded with parameters empirically derived through statistical analysis and optimizations with the SKJ model, see *put_larvae_at_obs* in *SeapodymCoupled_EarlyLife.cpp*. 

2) Small change in Dinf (maximal theoretical diffusion) definition - instead of using linear relationship with size given by line $y=Dspeed(1.25-.25 L(age))$ for the smallest to largest tuna, the allometric function analogous to the one used for advection rates is now used. For the moment, there still remain hard-coded parameters, so further revisions will be necessary (at least taking out parameters to the parfile), but the results of optimizations has shown some improvements in estimation of habitats yielding lower diffusion rates for largest tunas. See *precaldia_comp* function in *caldia.cpp*.

3) Fixed a minor bug in the Taglike function, due to which the first group of tags was never used in the likelihood.

4) Taylor test was added and essentially replaces the Autodif's derivative check. Current implementation has hard-coded number of 16 points for finite difference calculations, although it can be envisaged to enter user specified number and pass it as an argument to the function, see *Taylor_derivative_test* in hessian.cpp.

5) Small correction of degraded cells in fisheries data for a closer match with the model grid, as well restriction of the LF data sample sized at smaller smax, see *Readwrite_fisheries.cpp*. 

6) Fixes in command line options for *seapodym_habitats* and *seapodym_densities*.

7) Small bug fix in *Hessian_comp* function, see hessian.cpp.

8) Removed all forcing data changes in the code.

9) Some code cleaning.


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

  to see the running options of this SEAPODYM sub-model. 

  To compile all model applications at once, run the following batch file inside folder containing Makefiles

    [~/seapodym/]$ . all.bat

  For convenience, the path to the _~/seapodym/bin/_ directory can be added to variable PATH in the user's _~/.bashrc_. 


## Running SEAPODYM models
    
  Example #1. **Habitat**

  Go to the *example-configs* directory and test that the binaries execute nominally and the models generate expected outputs on your computer by executing a small habitat model example.
    
    [~]$ cd ~/seapodym/examples-configs/habitat    

  *Simulation run*

  To execute habitat model in simulation mode type the command with option *-s*
    
    [~/seapodym/examples-configs/habitat]$ seapodym\_habitats -s habitat.xml

  The simulation log can be compared with the one provided in the sim.out file. 
  
  *Optimization run*

  If option *-s* is omited, the application _seapodym\_habitats_ will start optimization. Note, due to much higher CPU and RAM demands compared to simulation mode, running optimization experiments is advised only after familializing yourself with the model and method basics, described in [Model Reference Manual](docs/manual/Seapodym_user_manual.pdf). 

  The optimization in _seapodym\_habitats_ aims to estimate habitat parameters by fitting the modelled habitat to an a priori known (e.g., previously estimated) habitat field. This can be convenient when only habitat parameters need to be (re-)estimated, e.g., in the twin experiments, or once the ocean forcing fields or other external parameters had been modified. For this example, "observations" were generated by running _seapodym\_habitats_ simulation with parameters in file _habitat\_input.xml_. The 3D (time and two-dimensional space) habitat field is written in binary file _msp\_spawning\_habitat\_input.dym_ (currently stored in directory _output_), which can be read with help of [dym](https://github.com/PacificCommunity/seapodym-codebase/tree/master/rtools/dym) R library. To test this application in optimization mode, run it starting with non-optimal parameters

    [~/seapodym/examples-configs/habitat]$ seapodym\_habitats habitat.xml

  The output should be the same as in the log file _opt.out_. The new parameters written in _newparfile.xml_ once the optimization converged to a minimum should be very close to that in habitat\_input.xml.  


  Example #2. **Skipjack**

  This branch will include the configuration for revised skipjack model with JRA55 forcing. Once the forcings are uploaded to the OSF data repository, the parfile will be updated... 
  
  Once downloaded, unzip and place the forcing files into a local directory without modifying the folder structure.

  Now open the skipjack\_F0.xml parfile in the preferred text editor and replace the ${SEAPODYM\_HOME} by the full absolute path to the unzipped *data* folder.  
    
    [~]$ cd example-configs/skipjack/

  *Simulation run*
    
  To run this optimal model in simulation mode type the following command

    [~/seapodym/example-configs/skipjack/]$ seapodym -s skipjack_F0.xml
    
  The simulation log can be compared with the one provided in sim\_F0.out file. The binary outputs will be written in folder _output/output\_F0_.

  Note that running this simulation with fishing requires fisheries data, which are not public. To inquire for the access to the data, please email <ofpdatarequest@spc.int>. Your request will be examined based on the eligibility of your organisation and your project. Then, the simulation with fishing can be run with the XML configuration file skipjack.xml, and the screen log compared with provided sim.out file.

  *Optimization run*

  It is possible to run optimization without fisheries data using _seapodym\_densities_ application. Make sure it has been built. If not, run the batch file to compile all five applications (see *Building from source* instructions above). Similarly to _seapodym\_habitats_, this application uses model outputs as "observations". Using this example configuration, make a copy of _skipjack\_F0.xml_ parfile

    [~/seapodym/example-configs/skipjack/]$ cp skipjack_F0.xml initparfile.xml

  In _initparfile.xml_ reduce the simulation time period to five years by modifying the _year_ attribute in _save_last_date_ node from 2010 to 1983. First run simulation with this parfile to generate pseudo-observations

    [~/seapodym/example-configs/skipjack/]$ seapodym -s initparfile.xml

  Make the input file with pseudo-observations by copying the binary output file _skj\_totbm.dym_ 
    
    [~/seapodym/example-configs/skipjack/]$ cp output/output_F0/skj_totbm.dym output/output_F0/skj_density_input.dym

  Now modify the _initparfile.xml_ to have non-optimal parameters. For example, test the optimization with perturbed predation mortality slope coefficient _Mp\_mean\_exp_, adding 0.01 to it. Now you can run optimization

  [~/seapodym/example-configs/skipjack/]$ seapodym\_densities initparfile.xml

  and compare the screen outputs with _opt.out_, and parameter estimates with those in _skipjack\_F0.xml_ parfile, which were used to produce pseudo-observations.  
  
     
## Documentation
For further information on how to use SEAPODYM for simulation studies and model developments see the [Model Reference Manual](docs/manual/Seapodym_user_manual.pdf). Doxygen generated [Code Documentation](docs/code-dox/codedoc_seapodym.pdf) can be consulted for a quick overview of the C++ numerical model code. 

## License
SEAPODYM is an open source project. The model code, associated tools and documentations are provided under the general terms of the three-clause BSD [LICENSE](LICENSE.md).
