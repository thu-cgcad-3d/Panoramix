# vis module
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})

file (GLOB VIS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB VIS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file (GLOB VIS_RESOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.qrc)	
file (GLOB VIS_UIS ${CMAKE_CURRENT_SOURCE_DIR}/*.ui)	

set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)
set(CMAKE_AUTORCC on)

add_library (${VIS_STATIC_LIB} STATIC 
	${VIS_SOURCES} 
	${VIS_HEADERS} 
	${VIS_RESOURCES}
	${VIS_UIS}
)
target_link_libraries (
	${VIS_STATIC_LIB} 
	${CORE_STATIC_LIB} 
	${Fundamental_LIBS} 
	${Qt_LIBS}
)