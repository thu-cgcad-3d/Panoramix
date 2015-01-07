
find_package(CUDA REQUIRED)

file (GLOB GPUQUICKSHIFT_SOURCE *.h *.cpp *.cu)

include_directories (${CUDA_TOOLKIT_INCLUDE})
cuda_add_library ("gpuquickshift" STATIC ${GPUQUICKSHIFT_SOURCE})
TARGET_LINK_LIBRARIES("gpuquickshift" ${CUDA_LIBRARIES})

set_property(TARGET "gpuquickshift" PROPERTY FOLDER "${PROJECT_NAME}.3rdparty")
