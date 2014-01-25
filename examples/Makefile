include ../config.mk

CPPFLAGS = $(EIGEN_CPPFLAGS) $(OPT_CPPFLAGS) -I..
LDFLAGS  = $(EIGEN_LDFLAGS)

all: clark-y-section vawt
	
clark-y-section: clark-y-section.cpp ../vortexje/libvortexje.a
	$(CXX) -o clark-y-section clark-y-section.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a	
	
vawt: vawt.cpp ../vortexje/libvortexje.a
	$(CXX) -o vawt vawt.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a
	
clean:
	rm -f clark-y-section vawt
