add_subdirectory(airfoils)

set(SRCS
    ellipse-generator.cpp)
	
set(HDRS
    ellipse-generator.hpp)

add_library(shape-generators OBJECT ${SRCS})

install(FILES ${HDRS} DESTINATION include/vortexje/shape-generators)
