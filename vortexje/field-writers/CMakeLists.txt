set(SRCS
    vtk-field-writer.cpp)
	
set(HDRS
    vtk-field-writer.hpp)

add_library(field-writers OBJECT ${SRCS})

install(FILES ${HDRS} DESTINATION include/vortexje/field-writers)
