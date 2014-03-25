include_directories (${CMAKE_CURRENT_SOURCE_DIR})

# core MODULE
set (CORE_MODULE_NAME ${PROJECT_NAME}.Core)
set (CORE_STATIC_LIB ${CORE_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${CORE_MODULE_NAME}")
set (CORE_SHARED_LIB ${CORE_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${CORE_MODULE_NAME}")

add_subdirectory("core")



# vis MODULE
set (VIS_MODULE_NAME ${PROJECT_NAME}.Vis)
set (VIS_STATIC_LIB ${VIS_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${VIS_MODULE_NAME}")
set (VIS_SHARED_LIB ${VIS_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${VIS_MODULE_NAME}")

add_subdirectory("vis")
