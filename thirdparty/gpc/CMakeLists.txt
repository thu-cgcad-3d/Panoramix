
set (GPC_SOURCE gpc.h gpc.c)

source_group ("gpc" FILES ${GPC_SOURCE})
add_library ("gpc" STATIC ${GPC_SOURCE})

set_property(TARGET "gpc" PROPERTY FOLDER "${PROJECT_NAME}.3rdparty")