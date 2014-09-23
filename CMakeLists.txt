cmake_minimum_required (VERSION 2.8.12)

# project name
set (PROJECT_NAME Panoramix)
string (TOUPPER ${PROJECT_NAME} PROJECT_NAME_UPPER) # upper case name

# project version
set (PROJECT_VERSION_MAJOR 0)
set (PROJECT_VERSION_MINOR 1)
set (PROJECT_VERSION ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})

set (CMAKE_ALLOW_LOOSE_CONSTRUCTS true)


project (${PROJECT_NAME})
enable_testing ()

### env setting ###
# set cmake module paths
list (APPEND CMAKE_MODULE_PATH 
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules/ 
	${CMAKE_CURRENT_SOURCE_DIR}/cmake/scripts/
)


### compiler setting ###
if (NOT DEFINED MSVC)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
endif ()


# global definitions
include (DefineCompilerMacros)


### 3rd party libs setting ###
# get fundamental libs
include (LocateFundamentalLibs)


# get Qt libs
set (Qt_DIR "" CACHE PATH "Qt root directory here")
set (Qt_MODULES_REQUIRED Core Gui Widgets OpenGL)
include (LocateQtLibs)


# set path env
include (DefineAddToEnvPath)
ADDTO_ENVPATH (Fundamental_PATH)
ADDTO_ENVPATH (Qt_PATH)


# include third party codes
add_subdirectory ("thirdparty")
list (APPEND Fundamental_INCLUDES ${ThirdParty_INCLUDES})
list (APPEND Fundamental_LIBS ${ThirdParty_LIBS})


# the libs
add_subdirectory (src lib)


# executable
add_subdirectory (executable executable)

