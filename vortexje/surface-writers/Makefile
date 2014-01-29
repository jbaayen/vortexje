include ../../config.mk

HPPDEPS = gmsh-surface-writer.hpp vtk-surface-writer.hpp

all:

install:
	mkdir -p $(PREFIX)/include/vortexje/surface-writers
	cp $(HPPDEPS) $(PREFIX)/include/vortexje/surface-writers

clean:
	rm -f *.o
