include config.mk

all:
	cd vortexje; $(MAKE)
	cd test; $(MAKE)
	
install:
	cd vortexje; $(MAKE) install
	mkdir -p $(PREFIX)/lib/pkgconfig
	cp vortexje.pc $(PREFIX)/lib/pkgconfig

clean:
	cd vortexje; $(MAKE) clean
	cd test; $(MAKE) clean
