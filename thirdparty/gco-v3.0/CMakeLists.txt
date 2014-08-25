
file (GLOB GCO_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.h ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

add_library ("gco-v3.0" STATIC 
	${GCO_SOURCES}
)

set_property(TARGET "gco-v3.0" PROPERTY FOLDER "${PROJECT_NAME}.3rdparty")