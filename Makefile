all: curved-patch

TRANSFINITE=/home/salvi/project/transfinite
TRIANGLE=/home/salvi/project/cl-nurbs/tests/shewchuk

INCLUDES=-I/usr/include/eigen3 -I$(TRIANGLE) \
         -I$(TRANSFINITE)/src/geom -I$(TRANSFINITE)/src/transfinite
LDFLAGS=-L$(TRANSFINITE)/debug/geom -L$(TRANSFINITE)/debug/transfinite
LDLIBS=-lgeom -ltransfinite -lm -lstdc++

CXXFLAGS=-std=c++17 -g -Wall $(INCLUDES)

OBJECTS=curved-patch.o \
	curved-gc.o \
	curved-cb.o \
        lsq-plane.o \
	harmonic.o \
	constrained-harmonic.o \
	curved-domain.o

curved-patch: $(OBJECTS) $(TRIANGLE)/triangle.o

.PHONY: clean
clean:
	$(RM) curved-patch $(OBJECTS)
