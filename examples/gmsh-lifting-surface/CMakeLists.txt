add_executable(gmsh-lifting-surface gmsh-lifting-surface.cpp)
target_link_libraries(gmsh-lifting-surface vortexje)

configure_file(gmsh-lifting-surface.msh ${CMAKE_CURRENT_BINARY_DIR}/gmsh-lifting-surface.msh COPYONLY)
