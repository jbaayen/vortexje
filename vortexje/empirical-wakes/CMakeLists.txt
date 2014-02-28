set(SRCS
    ramasamy-leishman-wake.cpp)
	
set(HDRS
    ramasamy-leishman-wake.hpp)

add_library(empirical-wakes OBJECT ${SRCS})

install(FILES ${HDRS} DESTINATION include/vortexje/empirical-wakes)
