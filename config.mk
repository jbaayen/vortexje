PREFIX = /usr/local

EIGEN_CPPFLAGS = `pkg-config --cflags eigen3`
EIGEN_LDFLAGS  =

PYTHON_CPPFLAGS = `pkg-config --cflags python`
PYTHON_LDFLAGS  = `pkg-config --libs python`

OPT_CPPFLAGS = -O3 -march=core2

CXX = g++ -fPIC
