include_directories (${CMAKE_CURRENT_SOURCE_DIR})

# deriv MODULE
set (DERIV_MODULE_NAME ${PROJECT_NAME}.Deriv)
set (DERIV_STATIC_LIB ${DERIV_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${DERIV_MODULE_NAME}")
set (DERIV_SHARED_LIB ${DERIV_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${DERIV_MODULE_NAME}")

add_subdirectory("deriv")


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


# sandbox MODULE
set (SANDBOX_MODULE_NAME ${PROJECT_NAME}.SandBox)
set (SANDBOX_STATIC_LIB ${SANDBOX_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${SANDBOX_MODULE_NAME}")
set (SANDBOX_SHARED_LIB ${SANDBOX_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${SANDBOX_MODULE_NAME}")

add_subdirectory("sandbox")
