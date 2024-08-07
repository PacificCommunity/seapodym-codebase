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

1.  **Precompiled binaries**

Supported operating systems are

-   Ubuntu
-   Red Hat Linux

Note, the compiler GNU C++ should be installed.

Download the precompiled binary release corresponding to your Linux system and gcc versions from [Releases](https://github.com/PacificCommunity/seapodym-codebase/releases/tag/seapodym-4.1/) directory. Extract contents of the zip file to a folder *\~/seapodym/*.

In the **Terminal** window, go inside the extracted *bin* folder and type

[\~/seapodym/bin/]\$ seapodym -h

This command should display the usage instruction and running options of SEAPODYM.

2.  **Building from source**

Source code compilation requires the GNU C++ compiler installed on the computer.

In addition, the following two libraries should be installed:

-   libxml2
-   AUTODIF

The libxml2 library is used to read and write all application parameters in the XML file. If not installed already, you can install it with the command:

```         
[~]$ sudo apt-get install libxml2-dev
```

choosing the package libxml2-dev, which contains header files that are required to compile SEAPODYM applications.

The AUTODIF libraries provide an array language extension to C ++ and can be installed as a part of the ADModel Builder software. Visit the [ADMB project](https://github.com/admb-project/admb) page and follow the instructions for [quick installation](https://github.com/admb-project/admb/blob/main/docs/install/QuickStartUnix.md) or [building from source](https://github.com/admb-project/admb/blob/main/docs/install/BuildingSourceUnix.md) of the latest release of the ADMB software.

Once ADMB software is installed, in the .bashrc file declare environment variable pointing to it, as follows:

```         
export ADMB_HOME=/your-path-to-admb-folder/
```

In the case when shared (dynamic) libraries are to be used, add the following line in the .bashrc:

```         
export LD_LIBRARY_PATH=$ADMB_HOME/lib:$LD_LIBRARY_PATH
```

Finally, create the following soft link for the optimized AUTODIF library in the *lib* folder that was compiled with the gcc version installed on your machine. If the directory containing ADMB package is ~*/admb/, the static library files are located in*~ /admb/lib/. The one that is required by *seapodym* is called \**libadmbo*\*.a, for example *libadmbo-x86_64-linux-g++11.a*. To create a symbolic link do

```         
[~/] cd ~/admb/lib

[~/admb/lib/]$ln -s libadmbo-x86_64-linux-g++11.a libadmbo.a
```

Once the required libraries have been installed and configured, proceed to installing SEAPODYM. The source code can be downloaded either from the [latest release](https://github.com/PacificCommunity/seapodym-codebase/releases/seapodym-4.1/), or, alternatively, checking out the github repository:

```         
[~]$ git clone https://github.com/PacificCommunity/seapodym-codebase.git
```

If downloaded a release, unpack the archive. Go into the folder with the *src* directory and Makefiles, which is presumably named *\~/seapodym/*. Compile one of the SEAPODYM applications by typing the command

```         
[~/seapodym/]$ make
```

If the compilation is successful, it will create the binary file *seapodym* in directory *bin*. This application runs the model in a simulation mode only. Executing

```         
[~/seapodym/bin/]$ seapodym -h
```

should display the usage instruction and running options. To compile the sub-model of species habitats type

```         
[~/seapodym/]$ make -f Makefile.hab
```

The binary file *seapodym_habitats* should be created. Type

```         
[~/seapodym/bin/]$ seapodym_habitats -h
```

to see the running options of this SEAPODYM sub-model.

To compile all model applications at once, run the following batch file inside folder containing Makefiles

```         
[~/seapodym/]$ . all.bat
```

For convenience, the path to the *\~/seapodym/bin/* directory can be added to variable PATH in the user's *\~/.bashrc*.

## Running SEAPODYM models

Example #1. **Habitat**

Go to the *example-configs* directory and test that the binaries execute nominally and the models generate expected outputs on your computer by executing a small habitat model example.

```         
[~]$ cd ~/seapodym/examples-configs/habitat    
```

*Simulation run*

To execute habitat model in simulation mode type the command with option *-s*

```         
[~/seapodym/examples-configs/habitat]$ seapodym\_habitats -s habitat.xml
```

The simulation log can be compared with the one provided in the sim.out file.

*Optimization run*

If option *-s* is omited, the application *seapodym_habitats* will start optimization. Note, due to much higher CPU and RAM demands compared to simulation mode, running optimization experiments is advised only after familializing yourself with the model and method basics, described in [Model Reference Manual](docs/manual/Seapodym_user_manual.pdf).

The optimization in *seapodym_habitats* can be used to estimate habitat parameters in two ways:

-   Fitting the modelled habitat to an a priori known (e.g., previously estimated) habitat field. This can be convenient when only habitat parameters need to be (re-)estimated, e.g., in the twin experiments, or once the ocean forcing fields or other external parameters had been modified. For this example, "observations" were generated by running *seapodym_habitats* simulation with parameters in file *habitat_input.xml*. The 3D (time and two-dimensional space) habitat field is written in binary file *msp_spawning_habitat_input.dym* (currently stored in directory *output*), which can be read with help of [dym](https://github.com/PacificCommunity/seapodym-codebase/tree/master/rtools/dym) R library. To test this application in optimization mode, run it starting with non-optimal parameters

-   Fitting the modelled habitat to an external source of data, which can be seasonally aggregated or not, and/or of categorical nature (i.e. available as potentially uneven intervals of values. If so, the xml parameter file must be modified as follows:

    -   \<fit_spawning_habitat_raw\> must be set to "0".

    -   \<q_sp_larvae\> must be indicated. Model predictions will be scaled by this constant coefficient before comparison to habitat input field.

    -   \<strdir_larvae\> and \<file_larvae_data\> should indicate the path to the habitat input DYM file.

    -   \<larvae_input_categorical\> can be set to "1" if habitat input is available as a categorical variable. If so, \<nb_larvae_cat\>, \<larvae_density_bins\> and \<larvae_density_last_bin_width\> must be specified.

    -   \<larvae_likelihood\> must be set to "1".

    -   The type of likelihood function to be used should be indicated using the \<larvae_likelihood_type\> field (available options: Gaussian Kernel, Poisson, Truncated Poisson, Zero-inflated Poisson, Zero-inflated Negative Binomial). Depending on the likelihood function, other fields such as \<likelihood_larvae_sigma\>, \<likelihood_larvae_beta\>, \<likelihood_larvae_probzero\> can be required.

    -   Other fields can be indicated, such as \<fit_null_larvae\>, \<weight_null_larvae\>.

        [\~/seapodym/examples-configs/habitat]\$ seapodym_habitats habitat.xml

The output should be the same as in the log file *opt.out*. The new parameters written in *newparfile.xml* once the optimization converged to a minimum should be very close to that in habitat_input.xml.

Example #2. **Skipjack**

This is a more comprehensive full model example, the pre-configured model of skipjack tuna, based on the [parameter estimation approach with fisheries and tagging data](https://cdnsciencepub.com/doi/full/10.1139/cjfas-2018-0470). To run this model, the corresponding forcing directory needs to be downloaded from the public [data repository](https://osf.io/h8u93) on the OSF platform.

Once downloaded, unzip and place the forcing files into a local directory without modifying the folder structure.

Now open the skipjack_F0.xml parfile in the preferred text editor and replace the \${SEAPODYM_HOME} by the full absolute path to the unzipped *data* folder.

```         
[~]$ cd example-configs/skipjack/
```

*Simulation run*

To run this optimal model in simulation mode type the following command

```         
[~/seapodym/example-configs/skipjack/]$ seapodym -s skipjack_F0.xml
```

The simulation log can be compared with the one provided in sim_F0.out file. The binary outputs will be written in folder *output/output_F0*.

Note that running this simulation with fishing requires fisheries data, which are not public. To inquire for the access to the data, please email [ofpdatarequest\@spc.int](mailto:ofpdatarequest@spc.int){.email}. Your request will be examined based on the eligibility of your organisation and your project. Then, the simulation with fishing can be run with the XML configuration file skipjack.xml, and the screen log compared with provided sim.out file.

*Optimization run*

It is possible to run optimization without fisheries data using *seapodym_densities* application. Make sure it has been built. If not, run the batch file to compile all five applications (see *Building from source* instructions above). Similarly to *seapodym_habitats*, this application uses model outputs as "observations". Using this example configuration, make a copy of *skipjack_F0.xml* parfile

```         
[~/seapodym/example-configs/skipjack/]$ cp skipjack_F0.xml initparfile.xml
```

In *initparfile.xml* reduce the simulation time period to five years by modifying the *year* attribute in *save_last_date* node from 2010 to 1983. First run simulation with this parfile to generate pseudo-observations

```         
[~/seapodym/example-configs/skipjack/]$ seapodym -s initparfile.xml
```

Make the input file with pseudo-observations by copying the binary output file *skj_totbm.dym*

```         
[~/seapodym/example-configs/skipjack/]$ cp output/output_F0/skj_totbm.dym output/output_F0/skj_density_input.dym
```

Now modify the *initparfile.xml* to have non-optimal parameters. For example, test the optimization with perturbed predation mortality slope coefficient *Mp_mean_exp*, adding 0.01 to it. Now you can run optimization

[\~/seapodym/example-configs/skipjack/]\$ seapodym_densities initparfile.xml

and compare the screen outputs with *opt.out*, and parameter estimates with those in *skipjack_F0.xml* parfile, which were used to produce pseudo-observations.

Example #3. **Albacore**

Similarly to skipjack example, this is the pre-configured example to run the albacore tuna model, [parameterized and validated with fisheries data](https://www.sciencedirect.com/science/article/pii/S0967064519301511). The configuration file will allow running the population dynamics model of the south Pacific stock. To run this model, download the forcing directory from the [OSF data repository](https://osf.io/j53hc).

Unzip and place the forcing files into a local directory without modifying the folder structure.

Open the albacore_F0.xml parfile and modify the \${SEAPODYM_HOME}, so to provide the full absolute path to the unzipped *data* folder.

```         
[~]$ cd example-configs/albacore/
```

and run

```         
[~/seapodym/example-configs/albacore/]$ seapodym -s albacore_F0.xml
```

The simulation log can be compared with the one provided in sim_F0.out file. Note that running this simulation with fishing requires fisheries data, which are not public. To inquire for the access to the data, please email [ofpdatarequest\@spc.int](mailto:ofpdatarequest@spc.int){.email}. Your request will be examined based on the eligibility of your organisation and your project. Then, the simulation with fishing can be run with the XML configuration file albacore.xml, and the screen log compared with provided sim.out file.

## Documentation

For further information on how to use SEAPODYM for simulation studies and model developments see the [Model Reference Manual](docs/manual/Seapodym_user_manual.pdf). Doxygen generated [Code Documentation](docs/code-dox/codedoc_seapodym.pdf) can be consulted for a quick overview of the C++ numerical model code.

## License

SEAPODYM is an open source project. The model code, associated tools and documentations are provided under the general terms of the three-clause BSD [LICENSE](LICENSE.md).
