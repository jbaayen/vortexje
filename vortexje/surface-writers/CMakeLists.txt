set(SRCS
    gmsh-surface-writer.cpp
    vtk-surface-writer.cpp)
	
set(HDRS
    gmsh-surface-writer.hpp
    vtk-surface-writer.hpp)

add_library(surface-writers OBJECT ${SRCS})

install(FILES ${HDRS} DESTINATION include/vortexje/surface-writers)
