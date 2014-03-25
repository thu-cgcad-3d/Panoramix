# vis module

set (CMAKE_AUTOMOC on)
add_definitions (${Qt_DEFINITIONS})
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})

file (GLOB VIS_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB VIS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file (GLOB VIS_UIS ${CMAKE_CURRENT_SOURCE_DIR}/*.ui)
file (GLOB VIS_RESOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.qrc)	

source_group ("vis" FILES ${VIS_HEADERS} ${VIS_SOURCES})
qt5_add_resources (RESOURCE_ADDED ${VIS_RESOURCES})

add_library (${VIS_STATIC_LIB} STATIC ${VIS_SOURCES} ${VIS_HEADERS} ${RESOURCE_ADDED})
target_link_libraries (${VIS_STATIC_LIB} ${CORE_STATIC_LIB} ${Fundamental_LIBS} ${Qt_LIBS})

qt5_use_modules(${VIS_STATIC_LIB} ${Qt_MODULES_REQUIRED})

