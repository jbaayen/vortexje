set(SRCS
    naca4-airfoil-generator.cpp)
	
set(HDRS
    naca4-airfoil-generator.hpp)

add_library(airfoils OBJECT ${SRCS})

install(FILES ${HDRS} DESTINATION include/vortexje/shape-generators/airfoils)
