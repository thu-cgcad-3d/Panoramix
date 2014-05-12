# deriv module

include_directories (${Fundamental_INCLUDES})

file (GLOB DERIV_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB DERIV_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

source_group ("deriv" FILES ${DERIV_HEADERS} ${DERIV_SOURCES})

add_library (${DERIV_STATIC_LIB} STATIC ${DERIV_SOURCES} ${DERIV_HEADERS})
target_link_libraries (${DERIV_STATIC_LIB} 
	${CORE_STATIC_LIB} 
	${Fundamental_LIBS}
)


