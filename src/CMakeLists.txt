
set_property(GLOBAL PROPERTY USE_FOLDERS On)

include_directories (${Fundamental_INCLUDES} ${MATLAB_INCLUDES})
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})

if (MSVC)
	add_definitions ( "/DMATLAB_CODE_DIR=\"${MATLAB_CODE_DIR}\"")
else ()
	add_definitions ( "-DMATLAB_CODE_DIR=\"${MATLAB_CODE_DIR}\"")
endif ()


file (GLOB CORE_HEADERS ./core/*.hpp)
file (GLOB CORE_SOURCES ./core/*.c*)

source_group ("core" FILES ${CORE_HEADERS} ${CORE_SOURCES})


file (GLOB ML_HEADERS ./ml/*.hpp)
file (GLOB ML_SOURCES ./ml/*.c*)

source_group ("ml" FILES ${ML_HEADERS} ${ML_SOURCES})


file (GLOB VIS_HEADERS ./vis/*.hpp)
file (GLOB VIS_SOURCES ./vis/*.c*)
file (GLOB VIS_RESOURCES ./vis/*.qrc)
file (GLOB VIS_UIS ./vis/*.ui)

set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)
set(CMAKE_AUTORCC on)

source_group ("vis" FILES ${VIS_HEADERS} ${VIS_SOURCES} ${VIS_RESOURCES} ${VIS_UIS})


file (GLOB MISC_HEADERS ./misc/*.hpp)
file (GLOB MISC_SOURCES ./misc/*.c*)

source_group ("misc" FILES ${MISC_HEADERS} ${MISC_SOURCES})


panoramix_add_library ( ${PROJECT_NAME} 
	${CORE_SOURCES} 
	${CORE_HEADERS}
	${ML_HEADERS} 
	${ML_SOURCES}
	${VIS_SOURCES} 
	${VIS_HEADERS} 
	${VIS_RESOURCES}
	${VIS_UIS}
	${MISC_HEADERS}
	${MISC_SOURCES}
)
target_link_libraries ( ${PROJECT_NAME} 
	${Fundamental_LIBS} 
	${MATLAB_LIBS}
	${Qt_LIBS}
)

set_property(TARGET ${PROJECT_NAME} PROPERTY FOLDER "${PROJECT_NAME}.Libs")