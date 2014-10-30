add_executable(test-elliptic-planform test-elliptic-planform.cpp)
target_link_libraries(test-elliptic-planform vortexje)

add_test(elliptic-planform test-elliptic-planform)

