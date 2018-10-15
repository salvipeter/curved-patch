all: curved-patch

TRANSFINITE=/home/salvi/project/transfinite
TRIANGLE=/home/salvi/project/cl-nurbs/tests/shewchuk

INCLUDES=-I/usr/include/eigen3 -I$(TRANSFINITE)/src/geom -I$(TRIANGLE)
LDFLAGS=-L$(TRANSFINITE)/debug/geom
LDLIBS=-lgeom -lm -lstdc++

CXXFLAGS=-std=c++17 -g -Wall $(INCLUDES)

OBJECTS=curved-patch.o \
        lsq-plane.o \
	harmonic.o

curved-patch: $(OBJECTS) $(TRIANGLE)/triangle.o

.PHONY: clean
clean:
	$(RM) curved-patch $(OBJECTS)
