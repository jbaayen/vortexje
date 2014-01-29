configure_file(sphere.msh ${CMAKE_CURRENT_BINARY_DIR}/sphere.msh COPYONLY)

add_executable(test-sphere test-sphere.cpp)
target_link_libraries(test-sphere vortexje)

add_test(sphere test-sphere)

configure_file(test-naca0012-reference.txt ${CMAKE_CURRENT_BINARY_DIR}/test-naca0012-reference.txt COPYONLY)

add_executable(test-naca0012 test-naca0012.cpp)
target_link_libraries(test-naca0012 vortexje)

add_test(naca0012 test-naca0012)
