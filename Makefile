NAME := bin/WebPlasmaDispersion.wt

# Defaults for dlg-hp.ornl.gov

CC := gcc
CPP := g++

MODULES := src include

CFLAGS := 
CPPFLAGS := -g -pg 

INCLUDEFLAGS := 
LIBS := 
LFLAGS := 

# Wt
WTDIR:=/home/dg6/code/wt/wt-3.2-gnu-4.6
INCLUDEFLAGS+= -I${WTDIR}/include 
LIBS+= -L${WTDIR}/lib -lwt -lwthttp -lboost_signals-mt -lboost_filesystem-mt -lboost_system-mt 

# O2csl
O2SCLDIR:=/home/dg6/code/o2scl
INCLUDEFLAGS+=-I${O2SCLDIR}/include
LIBS+= -L${O2SCLDIR}/lib -lo2scl

# Armadillo

LIBS+= -L/usr/lib64 -L/usr/lib64/atlas -lblas -lgslcblas -llapack -lclapack

# Matpack
MATPACKDIR:=/home/dg6/code/matpack/matpack
INCLUDEFLAGS+= -I${MATPACKDIR}/include
LIBS+= ${MATPACKDIR}/matpack.a

LINK := $(CPP) ${CPPFLAGS}

# You shouldn't have to go below here
#
# DLG: 	Added the -x c to force c file type so that 
# 		the .cu files will work too :)

DIRNAME = `dirname $1`
MAKEDEPS = $(CC) -MM -MG $2 -x c $3 | sed -e "s@^\(.*\)\.o:@.dep/$1/\1.d obj/$1/\1.o:@"

.PHONY : all

all : $(NAME)

# look for include files in each of the modules
INCLUDEFLAGS += $(patsubst %, -I%, $(MODULES))

CFLAGS += $(INCLUDEFLAGS)
CPPFLAGS += $(INCLUDEFLAGS) 

# determine the object files
SRCTYPES := c cpp 
OBJ := $(foreach srctype, $(SRCTYPES), $(patsubst %.$(srctype), obj/%.o, $(wildcard $(patsubst %, %/*.$(srctype), $(MODULES)))))

# link the program
$(NAME) : $(OBJ)
	$(LINK) $(LFLAGS) -o $@ $(OBJ) $(LIBS)
	@echo 'Run webApp using ...'
	@echo 'WT_TMP_DIR=/home/dg6/code/web-plasma-dispersion/tmp bin/WebPlasmaDispersion.wt --docroot ./ --http-address 0.0.0.0 --http-port 8080 -c ./wt_config.xml'

# calculate include dependencies
.dep/%.d : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CPPFLAGS), $<) > $@

obj/%.o : %.cpp
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CPP) $(CPPFLAGS) -c -o $@ $<

.dep/%.d : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.d$$||'`
	$(call MAKEDEPS,$(call DIRNAME, $<), $(CFLAGS), $<) > $@

obj/%.o : %.c
	@mkdir -p `echo '$@' | sed -e 's|/[^/]*.o$$||'`
	$(CC) $(CFLAGS) -c -o $@ $<

# include the C include dependencies
DEP := $(patsubst obj/%.o, .dep/%.d, $(OBJ))

ifneq ($(MAKECMDGOALS),clean)
include $(DEP)
endif

clean:
	-@rm $(NAME) $(OBJ) $(DEP) .dep/src/*

allclean: 
	-@rm $(NAME) $(OBJ) $(DEP) .dep/src/* webFace.wt src_webFace/*.o

