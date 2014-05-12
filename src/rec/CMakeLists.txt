# reconstruction module

include_directories (${Fundamental_INCLUDES})

file (GLOB REC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB REC_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

source_group ("rec" FILES ${REC_HEADERS} ${REC_SOURCES})

add_library (${REC_STATIC_LIB} STATIC ${REC_SOURCES} ${REC_HEADERS})
target_link_libraries (${REC_STATIC_LIB} 
	${CORE_STATIC_LIB} 
	${DERIV_STATIC_LIB}
	${VIS_STATIC_LIB} 
	${Fundamental_LIBS}
)