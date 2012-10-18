include ../config.mk

CPPDEPS = mesh.cpp solver.cpp wing.cpp wake.cpp collection.cpp wing-builder.cpp

HPPDEPS = mesh.hpp solver.hpp parameters.hpp wing.hpp wake.hpp collection.hpp wing-builder.hpp

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
	cp libvortexje.so.1.0.0 $(PREFIX)/lib/libvortexje.so.1.0.0
	ln -sf $(PREFIX)/lib/libvortexje.so.1.0.0 $(PREFIX)/lib/libvortexje.so.1
	ln -sf $(PREFIX)/lib/libvortexje.so.1.0.0 $(PREFIX)/lib/libvortexje.so
	mkdir -p $(PREFIX)/include/vortexje
	cp $(HPPDEPS) $(PREFIX)/include/vortexje

clean:
	rm -f *.o *.a *.so *.so.1 *.so.1.0.0
