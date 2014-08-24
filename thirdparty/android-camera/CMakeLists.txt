
file (GLOB MOSAIC_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/feature_mos/src/mosaic/*.h)
file (GLOB MOSAIC_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/feature_mos/src/mosaic/*.cpp)

# file (GLOB MOSAIC_RENDERER_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/feature_mos/src/mosaic_renderer/*.h)
# file (GLOB MOSAIC_RENDERER_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/feature_mos/src/mosaic_renderer/*.cpp)

source_group ("feature_mos/mosaic" FILES ${MOSAIC_HEADERS} ${MOSAIC_SOURCES})
# source_group ("feature_mos/mosaic_renderer" FILES ${MOSAIC_RENDERER_HEADERS} ${MOSAIC_RENDERER_SOURCES})

file (GLOB STAB_DB_VLVM_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/feature_stab/db_vlvm/*.h)
file (GLOB STAB_DB_VLVM_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/feature_stab/db_vlvm/*.cpp)

source_group ("feature_stab/db_vlvm" FILES ${STAB_DB_VLVM_HEADERS} ${STAB_DB_VLVM_SOURCES})

file (GLOB STAB_DB_REG_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/feature_stab/src/dbreg/*.h)
file (GLOB STAB_DB_REG_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/feature_stab/src/dbreg/*.cpp)

source_group ("feature_stab/db_reg" FILES ${STAB_DB_REG_HEADERS} ${STAB_DB_REG_SOURCES})

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/feature_mos/src" 
	"${CMAKE_CURRENT_SOURCE_DIR}/feature_stab/db_vlvm" 
	"${CMAKE_CURRENT_SOURCE_DIR}/feature_stab/src")

add_library ("android-camera" STATIC 
	${MOSAIC_HEADERS} ${MOSAIC_SOURCES}
	# ${MOSAIC_RENDERER_HEADERS} ${MOSAIC_RENDERER_SOURCES}
	${STAB_DB_VLVM_HEADERS} ${STAB_DB_VLVM_SOURCES}
	${STAB_DB_REG_HEADERS} ${STAB_DB_REG_SOURCES}
)

set_property(TARGET "android-camera" PROPERTY FOLDER "${PROJECT_NAME}.3rdparty")