add_executable(clark-y-section clark-y-section.cpp)
target_link_libraries(clark-y-section vortexje)

add_executable(naca0012-section-oscillation naca0012-section-oscillation.cpp)
target_link_libraries(naca0012-section-oscillation vortexje)

add_executable(vawt vawt.cpp)
target_link_libraries(vawt vortexje)
