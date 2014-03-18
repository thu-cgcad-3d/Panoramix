set (CMAKE_ALLOW_LOOSE_CONSTRUCTS true)

cmake_minimum_required (VERSION 2.8)

set (PROJECT_NAME Panoramix)
project (${PROJECT_NAME})
enable_testing ()

set (PROJECT_LIB_STATIC_NAME ${PROJECT_NAME}Static)
set (PROJECT_LIB_SHARED_NAME ${PROJECT_NAME})

set (PROJECT_VERSION_MAJOR 0)
set (PROJECT_VERSION_MINOR 1)
set (PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

set (CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
if (NOT DEFINED MSVC)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif ()

message ( STATUS "PROJECT_SOURCE_DIR: " ${PROJECT_SOURCE_DIR} )
message ( STATUS "CMAKE_MODULE_PATH: " ${CMAKE_MODULE_PATH} )
message ( STATUS "CMAKE_SYSTEM_PROCESSOR: " ${CMAKE_SYSTEM_PROCESSOR} )


# 3rd libs
find_package (Eigen3 REQUIRED)
find_package (OpenCV REQUIRED)

set (THIRDPARTY_INCLUDES ${EIGEN3_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS} ${OPENGL_INCLUDE_DIR})
foreach (i ${THIRDPARTY_INCLUDES})
	message (STATUS "3rd party includes: " ${i})
endforeach ()
set (THIRDPARTY_LIBS ${OpenCV_LIBS} ${OPENGL_LIBRARIES})
foreach (i ${THIRDPARTY_LIBS})
	message (STATUS "3rd party libs: " ${i})
endforeach ()
# add path to opencv dlls on win
if (DEFINED _OpenCV_LIB_PATH)
	set (THIRDPARTY_PATH ${_OpenCV_LIB_PATH})
endif ()


# google test
set (GTEST_DIR "gtest-1.7.0")


# global includes
include_directories (${THIRDPARTY_INCLUDES} ${GTEST_DIR}/include ${PROJECT_BINARY_DIR}/include)

add_subdirectory (src lib)
add_subdirectory (${GTEST_DIR})
add_subdirectory (tests bin)

