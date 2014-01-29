add_subdirectory(rply)

set(SRCS
    gmsh-surface-loader.cpp
    ply-surface-loader.cpp)
	
set(HDRS
    gmsh-surface-loader.hpp
    ply-surface-loader.hpp)

add_library(surface-loaders STATIC ${SRCS})

target_link_libraries(surface-loaders rply)

install(FILES ${HDRS} DESTINATION include/vortexje/surface-loaders)
