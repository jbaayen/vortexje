add_executable(test-vortex-core test-vortex-core.cpp)
target_link_libraries(test-vortex-core vortexje)

add_test(vortex-core test-vortex-core)
