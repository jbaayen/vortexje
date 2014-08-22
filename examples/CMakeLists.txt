add_executable(clarky-section clarky-section.cpp)
target_link_libraries(clarky-section vortexje)

configure_file(clarky.dat ${CMAKE_CURRENT_BINARY_DIR}/clarky.dat COPYONLY)

add_executable(naca0012-section-oscillation naca0012-section-oscillation.cpp)
target_link_libraries(naca0012-section-oscillation vortexje)

add_executable(naca0025-rounded-tips naca0025-rounded-tips.cpp)
target_link_libraries(naca0025-rounded-tips vortexje)

add_executable(gmsh-lifting-surface gmsh-lifting-surface.cpp)
target_link_libraries(gmsh-lifting-surface vortexje)

configure_file(gmsh-lifting-surface.msh ${CMAKE_CURRENT_BINARY_DIR}/gmsh-lifting-surface.msh COPYONLY)

add_executable(vawt vawt.cpp)
target_link_libraries(vawt vortexje)
