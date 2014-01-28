include ../config.mk

CPPDEPS = \
	surface.cpp \
	solver.cpp \
	parameters.cpp \
	lifting-surface.cpp \
	wake.cpp \
	body.cpp \
	surface-builder.cpp \
	lifting-surface-builder.cpp \
	shape-generators/ellipse-generator.cpp \
	shape-generators/airfoils/naca4-generator.cpp \
	shape-generators/airfoils/clark-y-generator.cpp

HPPDEPS = \
	surface.hpp \
	solver.hpp \
	parameters.hpp \
	lifting-surface.hpp \
	wake.hpp body.hpp \
	surface-builder.hpp \
	lifting-surface-builder.hpp

OBJS = $(CPPDEPS:.cpp=.o)

CPPFLAGS = $(EIGEN_CPPFLAGS) $(OPT_CPPFLAGS) -I..
LDFLAGS  = $(EIGEN_LDFLAGS)

all: libvortexje.a

libvortexje.a: $(OBJS) $(HPPDEPS)
	rm -f libvortexje.a
	ar rcs libvortexje.a $(OBJS)

libvortexje.so.1.0.0: $(OBJS) $(HPPDEPS)
	gcc -shared -Wl,-soname,libvortexje.so.1 -o libvortexje.so.1.0.0 $(OBJS) $(LDFLAGS)

install: libvortexje.so.1.0.0
	cd shape-generators; $(MAKE) install
	cp libvortexje.so.1.0.0 $(PREFIX)/lib/libvortexje.so.1.0.0
	ln -sf $(PREFIX)/lib/libvortexje.so.1.0.0 $(PREFIX)/lib/libvortexje.so.1
	ln -sf $(PREFIX)/lib/libvortexje.so.1.0.0 $(PREFIX)/lib/libvortexje.so
	mkdir -p $(PREFIX)/include/vortexje
	cp $(HPPDEPS) $(PREFIX)/include/vortexje

clean:
	cd shape-generators; $(MAKE) clean
	rm -f *.o *.a *.so *.so.1 *.so.1.0.0
