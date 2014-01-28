include ../../config.mk

HPPDEPS = ellipse.hpp

all:
	cd airfoils; $(MAKE)

install:
	cd airfoils; $(MAKE) install
	mkdir -p $(PREFIX)/include/vortexje/shapes
	cp $(HPPDEPS) $(PREFIX)/include/vortexje/shapes

clean:
	cd airfoils; $(MAKE) clean
	rm -f *.o
