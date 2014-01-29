add_subdirectory(surface-loaders)
add_subdirectory(surface-writers)
add_subdirectory(field-writers)
add_subdirectory(shape-generators)

set(SRCS
    surface.cpp 
	solver.cpp 
	parameters.cpp 
	lifting-surface.cpp 
	wake.cpp 
	body.cpp 
	surface-builder.cpp 
	lifting-surface-builder.cpp 
	surface-writer.cpp 
	field-writer.cpp)
	
set(HDRS
    surface.hpp 
	solver.hpp 
	parameters.hpp 
	lifting-surface.hpp 
	wake.hpp 
	body.hpp 
	surface-builder.hpp 
	lifting-surface-builder.hpp 
	surface-loader.hpp
	surface-writer.hpp 
	field-writer.hpp)

add_library(vortexje SHARED ${SRCS})

target_link_libraries(vortexje surface-loaders surface-writers field-writers shape-generators)

install (TARGETS vortexje DESTINATION lib)
install (FILES ${HDRS} DESTINATION include/vortexje)
