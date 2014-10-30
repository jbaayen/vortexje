add_executable(hawt hawt.cpp)
target_link_libraries(hawt vortexje)

configure_file(hawt-airfoil.dat ${CMAKE_CURRENT_BINARY_DIR}/hawt-airfoil.dat COPYONLY)
configure_file(hawt-blade.dat ${CMAKE_CURRENT_BINARY_DIR}/hawt-blade.dat COPYONLY)
