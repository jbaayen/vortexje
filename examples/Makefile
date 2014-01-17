include ../config.mk

CPPFLAGS = $(EIGEN_CPPFLAGS) $(OPT_CPPFLAGS) -I..
LDFLAGS  = $(EIGEN_LDFLAGS)

all: wing wing-builder vawt

wing: wing.cpp ../vortexje/libvortexje.a
	$(CXX) -o wing wing.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a
	
wing-builder: wing-builder.cpp ../vortexje/libvortexje.a
	$(CXX) -o wing-builder wing-builder.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a	
	
vawt: vawt.cpp ../vortexje/libvortexje.a
	$(CXX) -o vawt vawt.cpp $(CPPFLAGS) $(LDFLAGS) ../vortexje/libvortexje.a
	
clean:
	rm -f wing vawt wing-builder
