all:
	cd vortexje; $(MAKE)
	cd python; $(MAKE)
	cd test; $(MAKE)
	
install:
	cd vortexje; $(MAKE) install
	cd python; $(MAKE) install
	mkdir -p /usr/local/lib/pkgconfig
	cp vortexje.pc /usr/local/lib/pkgconfig

clean:
	cd vortexje; $(MAKE) clean
	cd test; $(MAKE) clean
