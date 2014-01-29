set(SRCS
    naca4-generator.cpp
    clark-y-generator.cpp)
	
set(HDRS
    naca4-generator.hpp
    clark-y-generator.hpp)

add_library(airfoils STATIC ${SRCS})

install(FILES ${HDRS} DESTINATION include/vortexje/shape-generators/airfoils)
