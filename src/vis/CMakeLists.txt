# vis module

add_definitions (${OGRE_DEFINITIONS})

include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES} ${OGRE_INCLUDES})

file (GLOB VIS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB VIS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
#file (GLOB VIS_UIS ${CMAKE_CURRENT_SOURCE_DIR}/*.ui)
file (GLOB VIS_RESOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.qrc)	


add_library (${VIS_STATIC_LIB} STATIC 
	${VIS_SOURCES} 
	${VIS_HEADERS} 
)
target_link_libraries (
	${VIS_STATIC_LIB} 
	${CORE_STATIC_LIB} 
	${Fundamental_LIBS} 
	${Qt_LIBS}
	${OGRE_LIBS}
)




