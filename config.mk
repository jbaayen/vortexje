PREFIX = /usr/local

EIGEN_CPPFLAGS = `pkg-config --cflags eigen3`
EIGEN_LDFLAGS  =

OPT_CPPFLAGS = -O3 -march=core2

CXX = g++ -fPIC
