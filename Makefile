HEADERS= \
XMLDocument2.h \
calpop.h \
Map.h \
Matrices.h \
Numfunc.h \
ReadWrite.h \
Param.h \
Date.h \
SaveTimeArea.h \
SeapodymCoupled.h \
SimtunaFunc.h \
mytypes.h \
SeapodymDocConsole.h \
ctrace.h \
Utilities.h \
VarMatrices.h \
VarParamCoupled.h \
VarSimtunaFunc.h \
NishikawaLike.h

SRCS= \
main_simulation.cpp \
ad_buffers.cpp \
XMLDocument2.cpp \
VarParamCoupled.cpp \
seapodym_coupled.cpp \
Map.cpp \
Matrices.cpp \
Numfunc.cpp \
ReadWrite_TXT.cpp \
ReadWrite_DYM.cpp \
ReadWrite_fisheries.cpp \
Param.cpp \
Date.cpp \
SaveTimeArea.cpp \
SimtunaFunc.cpp \
VarParamCoupled_xinit.cpp \
VarParamCoupled_reset.cpp \
SeapodymCoupled_Funcs.cpp \
SeapodymCoupled_Forage.cpp \
SeapodymCoupled_EditRunCoupled.cpp \
SeapodymCoupled_OnRunCoupled.cpp \
SeapodymCoupled_OnRunFirstStep.cpp \
SeapodymCoupled_OnReadForcing.cpp \
SeapodymCoupled_OnWriteOutput.cpp \
SeapodymCoupled_ReadTags.cpp \
SeapodymDocConsole_UpdateDisplay.cpp \
spawning_habitat.cpp \
juvenile_habitat.cpp \
mortality_sp.cpp \
caldia.cpp \
tridag_bet.cpp \
calrec_adre.cpp \
precalrec_juv.cpp \
calrec_precalrec.cpp \
predicted_catch.cpp \
predicted_catch_without_effort.cpp \
total_exploited_biomass.cpp \
total_obs_catch_age.cpp \
spawning.cpp \
accessibility.cpp \
feeding_habitat.cpp \
seasonal_switch.cpp \
total_mortality_comp.cpp \
food_requirement_index.cpp \
fd_spawning_habitat.cpp \
fd_juvenile_habitat.cpp \
fd_mortality_sp.cpp \
fd_caldia.cpp \
fd_tridag_bet.cpp \
fd_calrec_adre.cpp \
fd_survival.cpp \
fd_precalrec_juv.cpp \
fd_calrec_precalrec.cpp \
fd_predicted_catch.cpp \
fd_predicted_catch_without_effort.cpp \
fd_total_exploited_biomass.cpp \
fd_total_obs_catch_age.cpp \
fd_spawning.cpp \
fd_accessibility.cpp \
fd_feeding_habitat.cpp \
fd_seasonal_switch.cpp \
fd_total_pop.cpp \
fd_total_mortality_comp.cpp \
fd_food_requirement_index.cpp \
Calpop_caldia.cpp  \
Calpop_calrec.cpp  \
Calpop_InitCalPop.cpp  \
Calpop_precaldia.cpp  \
Calpop_precalrec.cpp  \
Calpop_tridag.cpp \
hessian.cpp \
like.cpp \
NishikawaLike.cpp

SRCPATH=DOM/src:src
INCPATH=-IDOM/src -Isrc
BINPATH=bin

##############################################################
OBJPATH=objs
ADMODEL=$(ADMB_HOME)

#DEBUG=-pg -g
DEBUG= -g
#production code
CFLAGS=-DTRUE=true -DFALSE=false -D __GNUDOS__ -Dlinux -O3 -DOPT_LIB -Wall -Wno-deprecated -I$(ADMODEL)/include -I/usr/include/libxml2 $(INCPATH)

###### with safe library ################################
#LFLAGS=-lm -L$(ADMODEL_HOME)/lib  -ladmb -lstdc++ -lxml2 

###### with optimized library ###########################
LFLAGS= -L$(ADMODEL)/lib  -ladmbo  -ldl -lstdc++ -lxml2 -lm

CC=gcc
LL=$(CC)


vpath %.cpp $(SRCPATH)
vpath %.h $(SRCPATH)

OBJECTS=$(SRCS:%.cpp=$(OBJPATH)/%.o)

##############################################################
export OBJECTS
export OBJPATH
export CFLAGS
export LFLAGS
export CC
export LL
##############################################################

all: init $(BINPATH)/seapodym

init:
	@test -d $(OBJPATH) || mkdir -v $(OBJPATH)
	@test -d $(BINPATH) || mkdir -v $(BINPATH)

test: init $(OBJECTS)
	make -f Makefile.test

docs: $(SRCS) $(HEADERS)
	@doxygen

$(BINPATH)/seapodym : $(OBJECTS)
	$(LL) -o$@ $(DEBUG) $^ $(LFLAGS)

$(OBJPATH)/%.o : %.cpp
	$(CC) -o$@ $(DEBUG) -c $(CFLAGS) $(filter %.cpp, $^)

$(OBJECTS) : $(HEADERS)

clean:
	@rm -vf $(OBJECTS) $(BINPATH)/seapodym
	@rm -vf $(OBJECTS) $(BINPATH)/seapodym.exe
	@rm -vf $(OBJECTS) $(BINPATH)/TestSeapodymCoupled
	@rm -vrf $(OBJPATH)
	#@rm -vrf docs/code-dox
	#@rm -vf gmon.out
