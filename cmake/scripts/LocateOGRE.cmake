# set OGRE

find_package (OGRE REQUIRED)
# OGRE_INCLUDE_DIR
# OGRE_LIBRARY_DBG
# OGRE_LIBRARY_REL
# OGRE_BINARY_DBG
# OGRE_BINARY_REL

set (OGRE_MODULES Overlay Paging RTShaderSystem Terrain Volume)

find_package (OIS)

list (GET OGRE_INCLUDE_DIR 0 OGRE_INCLUDES)
get_filename_component (OGRE_SDK_INCLUDE_DIR ${OGRE_INCLUDES} DIRECTORY) # get sdk include dir (parent folder of OGRE include)
set (OGRE_INCLUDES ${OGRE_SDK_INCLUDE_DIR} ${OIS_INCLUDE_DIR} ${OGRE_INCLUDES})
set (OGRE_LIBS 
	optimized ${OIS_LIBRARY_REL} debug ${OIS_LIBRARY_DBG}
	optimized ${OGRE_LIBRARY_REL} debug ${OGRE_LIBRARY_DBG})
get_filename_component (OGRE_SDK_DIR ${OGRE_SDK_INCLUDE_DIR} DIRECTORY)

foreach (i ${OGRE_MODULES})
	list (APPEND OGRE_INCLUDES ${OGRE_${i}_INCLUDE_DIR})
	list (APPEND OGRE_LIBS optimized ${OGRE_${i}_LIBRARY_REL} debug ${OGRE_${i}_LIBRARY_DBG})
endforeach ()

get_filename_component(OGRE_PATH_REL ${OGRE_BINARY_REL} DIRECTORY)
get_filename_component(OGRE_PATH_DBG ${OGRE_BINARY_DBG} DIRECTORY)
get_filename_component(OIS_PATH_REL ${OIS_BINARY_REL} DIRECTORY)
get_filename_component(OIS_PATH_DBG ${OIS_BINARY_DBG} DIRECTORY)

set(OGRE_PATH_REL ${OGRE_PATH_REL} ${OIS_PATH_REL})
set(OGRE_PATH_DBG ${OGRE_PATH_DBG} ${OIS_PATH_DBG})

message (STATUS "OGRE SDK dir: " ${OGRE_SDK_DIR})
foreach (i ${OGRE_INCLUDES})
	message (STATUS "OGRE includes: " ${i})
endforeach ()
foreach (i ${OGRE_LIBS})
	message (STATUS "OGRE libs: " ${i})
endforeach ()

set (OGRE_PLUGINS_FOLDER_DBG ${OGRE_SDK_DIR}/bin/debug)
set (OGRE_PLUGINS_FOLDER_REL ${OGRE_SDK_DIR}/bin/release)
set (OGRE_MEDIA_FOLDER ${OGRE_SDK_DIR}/media)

if (MSVC)
	set (OGRE_DEFINITIONS 
		"/DOGRE_PLUGINS_FOLDER_DBG=\"${OGRE_PLUGINS_FOLDER_DBG}\"" 
		"/DOGRE_PLUGINS_FOLDER_REL=\"${OGRE_PLUGINS_FOLDER_REL}\"" 
		"/DOGRE_MEDIA_FOLDER=\"${OGRE_MEDIA_FOLDER}\""
	)
else ()
	set (OGRE_DEFINITIONS 
		"-DOGRE_PLUGINS_FOLDER_DBG=\"${OGRE_PLUGINS_FOLDER_DBG}\"" 
		"-DOGRE_PLUGINS_FOLDER_REL=\"${OGRE_PLUGINS_FOLDER_REL}\"" 
		"-DOGRE_MEDIA_FOLDER=\"${OGRE_MEDIA_FOLDER}\""
	)
endif ()



