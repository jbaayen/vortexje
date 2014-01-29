set(SRCS
    naca4-airfoil-generator.cpp
    clark-y-airfoil-generator.cpp)
	
set(HDRS
    naca4-airfoil-generator.hpp
    clark-y-airfoil-generator.hpp)

add_library(airfoils OBJECT ${SRCS})

install(FILES ${HDRS} DESTINATION include/vortexje/shape-generators/airfoils)
