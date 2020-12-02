#
#  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
#  DATE WRITTEN:     Mar 2020
#  LAST MODIFIED:    Dec 2020
#  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
#  DESCRIPTION:      Compute Chebyshev coefficients from cspice states
#

CC = gcc
CPP = g++

CFLAGS = -I$(ODIR) -I$(IDIR)

ODIR = obj
CSPICE = ../../../../../Dropbox/CODE/cspice/lib
IDIR = include
SDIR = src

LIBS = -lm #-lnlopt

_DEPS1_inc = const.hpp c_functions.hpp chebyshev_coeff_gen.hpp chebyshev.hpp \
lsq_chebyshev_fit.hpp ephem_gen_chebyshev.hpp
DEPS1_inc = $(patsubst %,$(IDIR)/%,$(_DEPS1_inc))

_SPICE = cspice.a
SPICE_lib = $(patsubst %,$(CSPICE)/%,$(_SPICE))

_OBJ1 = main_cspice.o c_functions.o chebyshev_coeff_gen.o chebyshev.o \
lsq_chebyshev_fit.o ephem_gen_chebyshev.o
OBJ1 = $(patsubst %,$(ODIR)/%,$(_OBJ1))

$(ODIR)/%.o: $(SDIR)/%.cpp
	$(CPP) -c -o $@ $< $(CFLAGS)

$(ODIR)/%.o: %.c $(DEPS1_inc)
	$(CC) -c -o $@ $< $(CFLAGS)

main_cspice: $(OBJ1) $(SPICE_lib)
	$(CPP) -g -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(INCDIR)/*~
