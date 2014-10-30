add_executable(clarky-section clarky-section.cpp)
target_link_libraries(clarky-section vortexje)

configure_file(clarky.dat ${CMAKE_CURRENT_BINARY_DIR}/clarky.dat COPYONLY)
