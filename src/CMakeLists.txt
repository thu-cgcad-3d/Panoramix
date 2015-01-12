
set_property(GLOBAL PROPERTY USE_FOLDERS On)

include_directories (${Fundamental_INCLUDES} ${MATLAB_INCLUDES})
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})

if (MSVC)
	add_definitions ( "/DMATLAB_CODE_DIR=\"${MATLAB_CODE_DIR}\"")
else ()
	add_definitions ( "-DMATLAB_CODE_DIR=\"${MATLAB_CODE_DIR}\"")
endif ()


file (GLOB CORE_HEADERS ./core/*.hpp)
file (GLOB CORE_SOURCES ./core/*.cpp)

source_group ("core" FILES ${CORE_HEADERS} ${CORE_SOURCES})

file (GLOB VIS_HEADERS ./vis/*.hpp)
file (GLOB VIS_SOURCES ./vis/*.cpp)
file (GLOB VIS_RESOURCES ./vis/*.qrc)
file (GLOB VIS_UIS ./vis/*.ui)

source_group ("vis" FILES ${VIS_HEADERS} ${VIS_SOURCES} ${VIS_RESOURCES} ${VIS_UIS})

set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)
set(CMAKE_AUTORCC on)

add_library ( ${PROJECT_NAME} STATIC 
	${CORE_SOURCES} 
	${CORE_HEADERS}
	${VIS_SOURCES} 
	${VIS_HEADERS} 
	${VIS_RESOURCES}
	${VIS_UIS}
)
target_link_libraries ( ${PROJECT_NAME} 
	${Fundamental_LIBS} 
	${MATLAB_LIBS}
	${Qt_LIBS}
)

set_property(TARGET ${PROJECT_NAME} PROPERTY FOLDER "${PROJECT_NAME}.Libs")