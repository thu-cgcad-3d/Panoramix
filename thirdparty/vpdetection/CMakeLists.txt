# file (GLOB JLINKAGE_INCLUDE_FILES ./JLinkage/Include/bm/*.h ./JLinkage/Include/kdtree++/*.hpp ./JLinkage/Include/*.h)
file (GLOB JLINKAGE_SRC_FILES ./JLinkage/Src/*.h ./JLinkage/Src/*.cpp)
file (GLOB VPDETECTION_FILES ./*cpp ./*h)

include_directories (./JLinkage/Include)
add_library ("vpdetection" STATIC ${JLINKAGE_SRC_FILES} ${VPDETECTION_FILES})

set_property(TARGET "vpdetection" PROPERTY FOLDER "${PROJECT_NAME}.3rdparty")