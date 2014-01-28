include ../../../config.mk

HPPDEPS = clark-y-generator.hpp naca4-generator.hpp

all:

install:
	mkdir -p $(PREFIX)/include/vortexje/shapes/airfoils
	cp $(HPPDEPS) $(PREFIX)/include/vortexje/shapes/airfoils

clean:
	rm -f *.o
