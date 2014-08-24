# core module

include_directories (${Fundamental_INCLUDES})

file (GLOB CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB CORE_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

#source_group ("core" FILES ${CORE_HEADERS} ${CORE_SOURCES})

add_library (${CORE_STATIC_LIB} STATIC ${CORE_SOURCES} ${CORE_HEADERS})
target_link_libraries (${CORE_STATIC_LIB}
	${Fundamental_LIBS}
)

