add_subdirectory(boundary-layers)
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
	boundary-layer.hpp
	body.hpp 
	surface-builder.hpp 
	lifting-surface-builder.hpp 
	surface-loader.hpp
	surface-writer.hpp 
	field-writer.hpp)

add_library(vortexje SHARED ${SRCS}
    $<TARGET_OBJECTS:boundary-layers>
    $<TARGET_OBJECTS:surface-loaders>
    $<TARGET_OBJECTS:surface-writers>
    $<TARGET_OBJECTS:field-writers>
    $<TARGET_OBJECTS:shape-generators>
    $<TARGET_OBJECTS:airfoils>
    $<TARGET_OBJECTS:rply>)

install (TARGETS vortexje DESTINATION lib)
install (FILES ${HDRS} DESTINATION include/vortexje)
