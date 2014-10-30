configure_file(sphere.msh ${CMAKE_CURRENT_BINARY_DIR}/sphere.msh COPYONLY)

add_executable(test-sphere test-sphere.cpp)
target_link_libraries(test-sphere vortexje)

add_test(sphere test-sphere)
