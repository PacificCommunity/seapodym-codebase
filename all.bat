make
objcopy --strip-all bin/seapodym bin/seapodym

make -f Makefile.clt
objcopy --strip-all bin/seapodym_cltags bin/seapodym_cltags

make -f Makefile.flx
objcopy --strip-all bin/seapodym_fluxes bin/seapodym_fluxes

make -f Makefile.den
objcopy --strip-all bin/seapodym_densities bin/seapodym_densities

make -f Makefile.hab
objcopy --strip-all bin/seapodym_habitats bin/seapodym_habitats
