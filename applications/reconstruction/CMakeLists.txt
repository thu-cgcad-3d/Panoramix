set (APP_NAME "Reconstruction")

set (CMAKE_AUTOMOC on)
set (CMAKE_INCLUDE_CURRENT_DIR on)

add_definitions (${Qt_DEFINITIONS})
add_definitions (${OGRE_DEFINITIONS})

include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES} ${OGRE_INCLUDES})
include_directories (${CMAKE_CURRENT_BINARY_DIR}) 

file (GLOB APP_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB APP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
file (GLOB APP_RESOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.qrc)	
file (GLOB APP_UIS ${CMAKE_CURRENT_SOURCE_DIR}/*.ui)

#message (STATUS "added resources: " ${RESOURCE_ADDED})

qt5_wrap_ui (UI_ADDED ${APP_UIS})
#qt5_add_resources (RESOURCE_ADDED ${APP_RESOURCES})

set_source_files_properties(${RESOURCE_ADDED} ${UI_ADDED} PROPERTIES GENERATED on)


# define app data location
get_filename_component(DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/" REALPATH)
message (STATUS "app data directory: " ${DATA_DIR})
if (MSVC)
	add_definitions ( "/DPROJECT_DATA_DIR_STR=\"${DATA_DIR}\"" )
else ()
	add_definitions ( "-DPROJECT_DATA_DIR_STR=\"${DATA_DIR}\"" )
endif ()

if (MSVC)
	list (REMOVE_DUPLICATES OGRE_PATH_DBG)
	list (REMOVE_DUPLICATES OGRE_PATH_REL)
	set (OGRE_PATH ${OGRE_PATH_DBG} ${OGRE_PATH_REL})
	list (REMOVE_DUPLICATES OGRE_PATH)
	JOIN("${Fundamental_PATH};${Qt_PATH};${OGRE_PATH}" ";" LD_PATH)

	configure_file (${CMAKE_CURRENT_SOURCE_DIR}/setpath.bat.in ${CMAKE_CURRENT_BINARY_DIR}/setpath.bat)
	message (WARNING "CMake cannot configure Visual Studio to modify the environment path during program execution, "
		"therefore you have to do this manually: add 'PATH=\$(PATH);${LD_PATH};' to [Project Property]->[Debug]->[Environment] ")
endif ()


# configure files for OGRE
string(REGEX REPLACE "/" "\\\\" OGRE_PLUGINS_FOLDER_DBG ${OGRE_PLUGINS_FOLDER_DBG}) 
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/resources.cfg.in ${CMAKE_CURRENT_BINARY_DIR}/resources.cfg)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/plugins.cfg.in ${CMAKE_CURRENT_BINARY_DIR}/plugins.cfg)
configure_file (${CMAKE_CURRENT_SOURCE_DIR}/plugins_d.cfg.in ${CMAKE_CURRENT_BINARY_DIR}/plugins_d.cfg)


add_executable (${APP_NAME} ${APP_SOURCES} ${APP_HEADERS} ${UI_ADDED} ${RESOURCE_ADDED})
target_link_libraries (${APP_NAME} 
	${CORE_STATIC_LIB} ${DERIV_STATIC_LIB} ${VIS_STATIC_LIB} ${REC_STATIC_LIB} 
	${Fundamental_LIBS} ${Qt_LIBS} ${OGRE_LIBS})


set_property(GLOBAL PROPERTY USE_FOLDERS On)
set_property(TARGET ${APP_NAME} PROPERTY FOLDER "${PROJECT_NAME}.Applications")




