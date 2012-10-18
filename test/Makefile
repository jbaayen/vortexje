include ../config.mk

CPPFLAGS = $(EIGEN_CPPFLAGS) $(OPT_CPPFLAGS) -I..
LDFLAGS  = $(EIGEN_LDFLAGS)

all: test-wing test-wing-builder test-vawt

test-wing: test-wing.cpp ../vortexje/libvortexje.a
	$(CXX) -o test-wing test-wing.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a
	
test-wing-builder: test-wing-builder.cpp ../vortexje/libvortexje.a
	$(CXX) -o test-wing-builder test-wing-builder.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a	
	
test-vawt: test-vawt.cpp ../vortexje/libvortexje.a
	$(CXX) -o test-vawt test-vawt.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a
	
clean:
	rm -f test-wing test-vawt test-wing-builder
