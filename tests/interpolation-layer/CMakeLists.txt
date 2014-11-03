add_executable(test-interpolation-layer test-interpolation-layer.cpp)
target_link_libraries(test-interpolation-layer vortexje)

add_test(interpolation-layer test-interpolation-layer)
