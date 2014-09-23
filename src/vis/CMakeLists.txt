# vis module
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})

file (GLOB VIS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB VIS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file (GLOB VIS_RESOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.qrc)	


add_library (${VIS_STATIC_LIB} STATIC 
	${VIS_SOURCES} 
	${VIS_HEADERS} 
	${VIS_RESOURCES}
)
target_link_libraries (
	${VIS_STATIC_LIB} 
	${CORE_STATIC_LIB} 
	${Fundamental_LIBS} 
	${Qt_LIBS}
)





