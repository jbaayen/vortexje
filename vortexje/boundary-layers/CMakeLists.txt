set(SRCS
    dummy-boundary-layer.cpp)
	
set(HDRS
    dummy-boundary-layer.hpp)

add_library(boundary-layers OBJECT ${SRCS})

install(FILES ${HDRS} DESTINATION include/vortexje/boundary-layers)
