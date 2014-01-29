add_subdirectory(airfoils)

set(SRCS
    ellipse-generator.cpp)
	
set(HDRS
    ellipse-generator.hpp)

add_library(shape-generators STATIC ${SRCS})

target_link_libraries(shape-generators airfoils)

install(FILES ${HDRS} DESTINATION include/vortexje/shape-generators)
