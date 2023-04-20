# Installation

## 1. Preliminaries

A list of dependencies and their sub-dependencies is presented below :
- [GCC](https://gcc.gnu.org/)
- [G++](https://gcc.gnu.org/) (included in GCC - Gnu Compiler Collection)
- [CMake](https://cmake.org/)
- [Libxml2](https://gitlab.gnome.org/GNOME/libxml2)
- [ADMB](https://www.admb-project.org/)
    - [Flex](https://github.com/westes/flex)

> Note : Many compilation errors can be caused by using a too old version of these packages. It is therefore advisable to keep these packages as up-to-date as possible.
>
>If you are not able to update them (e.g. if you are not the system administrator), it is recommended to use the Anaconda environment (see also [A. With Anaconda](#a-with-anaconda)).

### A. With [Anaconda](https://docs.anaconda.com/)

First of all, Anaconda must be installed on your system. Select your system [here](https://docs.anaconda.com/anaconda/install/) and follow the instructions.

Once Anaconda is installed, you have to create an environment :

>**Flex** is required to compile **ADMB** from source code.
>
>Fairly recent versions of **Gcc**, **Gxx** and **CMake** are required to compile the **SEAPODYM** project.
```bash
conda create -n seapodym -c conda-forge gcc gxx cmake flex libxml2
```

After what you have to activate your new environment :
```bash
conda activate seapodym
```

### B. Without Anaconda

Install the needed package using your favorit package manager (apt, apt-get, snap etc...).
```bash
sudo apt update
sudo apt full-upgrade
sudo apt install gcc cmake flex libxml2 libxml2-dev
```

## 2. [ADMB](https://www.admb-project.org/)
The **ADMB** project supports the application of automatic differentiation (AD) for solutions to non-linear statistical modeling and optimization problems.

**Autodif** libraries provide an array language extension to C ++ enabling the automatic code differentiation ([*Autodif User’s Manual*](https://www.admb-project.org/docs/manuals/)).

If it is not installed already, you can download and build the library from [official releases](http://www.admb-project.org/news/) or install the latest Github release following these steps :

- Clone the [ADMB GitHub repository](https://github.com/admb-project/admb).
```bash
git clone https://github.com/admb-project/admb.git
```
- Compile the code using the makefile.
> `make` if you don't need the shared libraries. `make all` instead.
```bash
cd admb
make all
```


- *Optional* - Test the installation.
```bash
make --directory=examples all
```

- *Optional* - If you have the admin rights you can move the files to the `/usr/local/admb` directory.
```bash
sudo make install
```
**But** if you don't, the build path (`/admb-folder-path/build`) will be used in CMake to compile **SEAPODYM** project.

## 3. [SEAPODYM](https://github.com/Ash12H/seapodym-codebase)

**Seapodym** use the **CMake** interface to build the executable files. Follow these steps :

- Clone the Seapodym repository or download and extract the project from Github :
```bash
git clone https://github.com/Ash12H/seapodym-codebase.git
```
- Create a build folder to store all the CMake compilation files :
```bash
cd seapodym-codebase
mkdir build
cd build
```

- Now their is two possibilities. If you installed **ADMB** using the `sudo make install` command, you can simply execute cmake configuration that way :
```bash
cmake ..
```
**But** if you did not installed ADMB in standard directories, you must pass corresponding arguments :
- `-D CMAKE_PREFIX_PATH=/path_to_admb_folder/build/`

> Since ADMB creates libraries which name contains informations about compiler and system configuration, it is recommanded to rename both static (`.a`) and dynamic (`.so`) libraries this way :
> - `libadmb.a` and `libadmb.so`
> - `libadmbo.a` and `libadmbo.so`
> 
> Otherwise, you will need to specify the ADMB library you want to use to compile SEAPODYM with this argument : `-D ADMB_LIBRARY=/full_path_to_admb_library/admb_library_name.a` or `.so`

> **Warning** : The libraries we talked about in the previous section are not the ones that contain the keywords "contrib" or "shared".

- Once the configuration is done, you can produce the executable files :
```bash
make
```

Files are now available in the `seapodym-codebase/build/bin` directory.
