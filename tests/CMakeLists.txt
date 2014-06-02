configure_file(sphere.msh ${CMAKE_CURRENT_BINARY_DIR}/sphere.msh COPYONLY)

add_executable(test-sphere test-sphere.cpp)
target_link_libraries(test-sphere vortexje)

add_test(sphere test-sphere)

configure_file(naca0012-reference-data.txt ${CMAKE_CURRENT_BINARY_DIR}/naca0012-reference-data.txt COPYONLY)

add_executable(test-naca0012-airfoil test-naca0012-airfoil.cpp)
target_link_libraries(test-naca0012-airfoil vortexje)

add_test(naca0012 test-naca0012-airfoil)

add_executable(test-elliptic-planform test-elliptic-planform.cpp)
target_link_libraries(test-elliptic-planform vortexje)

add_test(elliptic-planform test-elliptic-planform)

