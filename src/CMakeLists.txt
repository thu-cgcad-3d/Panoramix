
set_property(GLOBAL PROPERTY USE_FOLDERS On)

include_directories (${Fundamental_INCLUDES} ${MATLAB_INCLUDES})
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})

if (MSVC)
	add_definitions ( "/DMATLAB_CODE_DIR=\"${MATLAB_CODE_DIR}\"")
	add_definitions ( "/DPROJECT_CACHE_DIR_STR=\"${CACHE_DIR}\"" )
else ()
	add_definitions ( "-DMATLAB_CODE_DIR=\"${MATLAB_CODE_DIR}\"")
	add_definitions ( "-DPROJECT_CACHE_DIR_STR=\"${CACHE_DIR}\"" )
endif ()



file (GLOB CORE_HEADERS ./core/*.hpp)
file (GLOB CORE_SOURCES ./core/*.c*)

source_group ("core" FILES ${CORE_HEADERS} ${CORE_SOURCES})


file (GLOB ML_HEADERS ./ml/*.hpp)
file (GLOB ML_SOURCES ./ml/*.c*)

source_group ("ml" FILES ${ML_HEADERS} ${ML_SOURCES})


file (GLOB GUI_HEADERS ./gui/*.hpp)
file (GLOB GUI_SOURCES ./gui/*.c*)
file (GLOB GUI_RESOURCES ./gui/*.qrc)
file (GLOB GUI_UIS ./gui/*.ui)

set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)
set(CMAKE_AUTORCC on)

source_group ("gui" FILES ${GUI_HEADERS} ${GUI_SOURCES} ${GUI_RESOURCES} ${GUI_UIS})


file (GLOB MISC_HEADERS ./misc/*.hpp)
file (GLOB MISC_SOURCES ./misc/*.c*)

source_group ("misc" FILES ${MISC_HEADERS} ${MISC_SOURCES})

file (GLOB EXPERIMENTAL_HEADERS ./experimental/*.hpp)
file (GLOB EXPERIMENTAL_SOURCES ./experimental/*.c*)

source_group ("experimental" FILES ${EXPERIMENTAL_HEADERS} ${EXPERIMENTAL_SOURCES})


panoramix_add_library ( ${PROJECT_NAME} 
	${CORE_SOURCES} 
	${CORE_HEADERS}
	${ML_HEADERS} 
	${ML_SOURCES}
	${GUI_SOURCES} 
	${GUI_HEADERS} 
	${GUI_RESOURCES}
	${GUI_UIS}
	${MISC_HEADERS}
	${MISC_SOURCES}
	${EXPERIMENTAL_HEADERS}
	${EXPERIMENTAL_SOURCES}
)
target_link_libraries ( ${PROJECT_NAME} 
	${Fundamental_LIBS} 
	${MATLAB_LIBS}
	${Qt_LIBS}
)

set_property(TARGET ${PROJECT_NAME} PROPERTY FOLDER "${PROJECT_NAME}.Libs")