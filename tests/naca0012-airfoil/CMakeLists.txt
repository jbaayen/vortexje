configure_file(naca0012-reference-data.txt ${CMAKE_CURRENT_BINARY_DIR}/naca0012-reference-data.txt COPYONLY)

add_executable(test-naca0012-airfoil test-naca0012-airfoil.cpp)
target_link_libraries(test-naca0012-airfoil vortexje)

add_test(naca0012 test-naca0012-airfoil)

