
set_property(GLOBAL PROPERTY USE_FOLDERS On)


# core MODULE
set (CORE_MODULE_NAME Core)
set (CORE_STATIC_LIB ${CORE_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${CORE_MODULE_NAME}")
set (CORE_SHARED_LIB ${CORE_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${CORE_MODULE_NAME}")

add_subdirectory("core")
set_property(TARGET ${CORE_STATIC_LIB} PROPERTY FOLDER "${PROJECT_NAME}.Libs")


# deriv MODULE
set (DERIV_MODULE_NAME Deriv)
set (DERIV_STATIC_LIB ${DERIV_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${DERIV_MODULE_NAME}")
set (DERIV_SHARED_LIB ${DERIV_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${DERIV_MODULE_NAME}")

add_subdirectory("deriv")
set_property(TARGET ${DERIV_STATIC_LIB} PROPERTY FOLDER "${PROJECT_NAME}.Libs")

# vis MODULE
set (VIS_MODULE_NAME Vis)
set (VIS_STATIC_LIB ${VIS_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${VIS_MODULE_NAME}")
set (VIS_SHARED_LIB ${VIS_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${VIS_MODULE_NAME}")

add_subdirectory("vis")
set_property(TARGET ${VIS_STATIC_LIB} PROPERTY FOLDER "${PROJECT_NAME}.Libs")


# rec MODULE
set (REC_MODULE_NAME Rec)
set (REC_STATIC_LIB ${REC_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${REC_MODULE_NAME}")
set (REC_SHARED_LIB ${REC_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${REC_MODULE_NAME}")

add_subdirectory("rec")
set_property(TARGET ${REC_STATIC_LIB} PROPERTY FOLDER "${PROJECT_NAME}.Libs")


# sandbox MODULE
set (SANDBOX_MODULE_NAME SandBox)
set (SANDBOX_STATIC_LIB ${SANDBOX_MODULE_NAME}.Static CACHE STRING 
	"static lib name of module ${SANDBOX_MODULE_NAME}")
set (SANDBOX_SHARED_LIB ${SANDBOX_MODULE_NAME}.Shared CACHE STRING
	"shared lib name of module ${SANDBOX_MODULE_NAME}")

add_subdirectory("sandbox")
set_property(TARGET ${SANDBOX_STATIC_LIB} PROPERTY FOLDER "${PROJECT_NAME}.Libs")
