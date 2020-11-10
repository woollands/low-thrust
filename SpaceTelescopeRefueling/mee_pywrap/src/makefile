#
#  AUTHORS:          Robyn Woollands (robyn.woollands@gmail.com)
#  DATE WRITTEN:     Oct 2019
#  LAST MODIFIED:    Oct 2019
#  AFFILIATION:      Jet Propulsion Laboratory, California Institute of Technology, Pasadena, CA
#  DESCRIPTION:      File to compile Picard-Chebyshev MEE code
#

# CC=clang
# CFLAGS=-Xpreprocessor -fopenmp -lomp
#
# obj= test_mee.o
#
# test: $(obj)
# 	 $(CC) ${CFLAGS} -o test $(obj)
#
# %.o:%.c
# 	 $(CC) ${CFLAGS} -c $<
#
# .PHONY: clean
#
# clean:
# 		rm *.o

IDIR = ../include
# CC = clang
CC = g++
# CC = gcc -std=c99 -g
# CC = gcc-9 -std=c99 -g -fopenmp -lomp
# CC = /usr/local/opt/llvm/bin/clang -fopenmp -L/usr/local/opt/llvm/lib gcc -std=c99 -g
# CC = /usr/local/opt/llvm/bin/clang -fopenmp gcc -std=c99 -g
# CC = clang -Xpreprocessor -fopenmp -lomp
# clang -Xpreprocessor -fopenmp -lomp test_mee.c -o test

CFLAGS = -I$(IDIR)
# CFLAGS = -I$(IDIR) -Xpreprocessor -fopenmp -lomp

ODIR = obj
LDIR = ../lib

LIBS = -lm

_DEPS1 = const.h c_functions.h classical2mee.h matrix_loader.h mee2rv.h ecef2eci.h eci2ecef.h EGM2008.h \
inertial2radial.h mee_state_costates.h
DEPS1 = $(patsubst %,$(IDIR)/%,$(_DEPS1))

_DEPS2 = clenshaw_curtis_ivpI.h lsq_chebyshev_fit.h chebyshev.h c_functions.h
DEPS2 = $(patsubst %,$(IDIR)/%,$(_DEPS2))

_OBJ1 = test_mee.o c_functions.o classical2mee.o matrix_loader.o mee2rv.o ecef2eci.o eci2ecef.o EGM2008.o \
inertial2radial.o mee_state_costates.o
OBJ1 = $(patsubst %,$(ODIR)/%,$(_OBJ1))

_OBJ2 = matrix_builder.o clenshaw_curtis_ivpI.o lsq_chebyshev_fit.o chebyshev.o c_functions.o clenshaw_curtis_ivpI.o
OBJ2 = $(patsubst %,$(ODIR)/%,$(_OBJ2))

_OBJ3 = mee_prop.o c_functions.o classical2mee.o matrix_loader.o mee2rv.o ecef2eci.o eci2ecef.o EGM2008.o \
inertial2radial.o mee_state_costates.o
OBJ3 = $(patsubst %,$(ODIR)/%,$(_OBJ3))

$(ODIR)/%.o: %.c $(DEPS1) $(DEPS2)
	$(CC) -c -o $@ $< $(CFLAGS)

test_mee: $(OBJ1)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LIBS)

matrix_builder: $(OBJ2)
	$(CC) -g -o $@ $^ $(CFLAGS) $(LIBS)

mee_prop: $(OBJ3)
	$(LINK.c) -fPIC -shared -o mee_prop.so $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.so $(ODIR)/*.o *~ core $(INCDIR)/*~

	# g++ -c -fPIC mee_prop.cpp matrix_loader.c mee2rv.c ecef2eci.c eci2ecef.c EGM2008.c inertial2radial.c c_functions.c
	# g++ mee_prop.o matrix_loader.o mee2rv.o ecef2eci.o eci2ecef.o EGM2008.o inertial2radial.o c_functions.o -shared -o mee_prop.so
