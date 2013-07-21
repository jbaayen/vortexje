include ../config.mk

CPPFLAGS = $(EIGEN_CPPFLAGS) $(OPT_CPPFLAGS) -I..
LDFLAGS  = $(EIGEN_LDFLAGS)

all: test-naca0012

test-naca0012: test-naca0012.cpp ../vortexje/libvortexje.a
	$(CXX) -o test-naca0012 test-naca0012.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a
	
clean:
	rm -f test-naca0012
