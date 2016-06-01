file (GLOB SRC "${CMAKE_CURRENT_SOURCE_DIR}/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/*c")
add_library ("wclique" STATIC ${SRC})
set_property(TARGET "wclique" PROPERTY FOLDER "${PROJECT_NAME}.3rdparty")