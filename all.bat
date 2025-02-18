make
objcopy --strip-all bin/seapodym bin/seapodym

make -f Makefile.clte
objcopy --strip-all bin/seapodym_clte bin/seapodym_clte

make -f Makefile.flx
objcopy --strip-all bin/seapodym_fluxes bin/seapodym_fluxes

make -f Makefile.den
objcopy --strip-all bin/seapodym_densities bin/seapodym_densities

make -f Makefile.hab
objcopy --strip-all bin/seapodym_habitats bin/seapodym_habitats
