
set_property(GLOBAL PROPERTY USE_FOLDERS On)

file (GLOB core_sources 
    "${CMAKE_CURRENT_SOURCE_DIR}/core/*.cpp" 
    "${CMAKE_CURRENT_SOURCE_DIR}/core/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/*.cu" 
    "${CMAKE_CURRENT_SOURCE_DIR}/core/*.cuh"
)
source_group ("core" FILES ${core_sources})

file (GLOB ml_sources 
    "${CMAKE_CURRENT_SOURCE_DIR}/ml/*.cpp" 
    "${CMAKE_CURRENT_SOURCE_DIR}/ml/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/ml/*.cu" 
    "${CMAKE_CURRENT_SOURCE_DIR}/ml/*.cuh"
)
source_group ("ml" FILES ${ml_sources})

file (GLOB gui_sources 
    "${CMAKE_CURRENT_SOURCE_DIR}/gui/*.cpp" 
    "${CMAKE_CURRENT_SOURCE_DIR}/gui/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/gui/*.cu" 
    "${CMAKE_CURRENT_SOURCE_DIR}/gui/*.cuh"
    "${CMAKE_CURRENT_SOURCE_DIR}/gui/*.qrc"
    "${CMAKE_CURRENT_SOURCE_DIR}/gui/*.ui"
)
set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)
set(CMAKE_AUTORCC on)
source_group ("gui" FILES ${gui_sources})

file (GLOB misc_sources 
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/*.cpp" 
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/*.cu" 
    "${CMAKE_CURRENT_SOURCE_DIR}/misc/*.cuh"
)
source_group ("misc" FILES ${misc_sources})

file (GLOB experimental_sources 
    "${CMAKE_CURRENT_SOURCE_DIR}/experimental/*.cpp" 
    "${CMAKE_CURRENT_SOURCE_DIR}/experimental/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/experimental/*.cu" 
    "${CMAKE_CURRENT_SOURCE_DIR}/experimental/*.cuh"
)
source_group ("experimental" FILES ${experimental_sources})



set (panoramix_sources 
	${core_sources} 
	${ml_sources} 
	${gui_sources}
	${misc_sources}
	${experimental_sources})
file (GLOB_RECURSE panoramix_test_sources "." *.test.cpp *.test.hpp *.test.cu *.test.cuh)
# remove test sources from the lib sources
list (REMOVE_ITEM panoramix_sources ${panoramix_test_sources})


if (MSVC)
    message (WARNING "CMake cannot configure Visual Studio to modify the environment path during program execution, "
        "therefore you have to do this manually: add 'PATH=\$(PATH);${DEPENDENCY_BINS_PATHS};' "
        "to [Project Property]->[Debug]->[Environment] ")
endif ()


# the lib project
message(STATUS "panoramix_sources:")
foreach(i ${panoramix_sources})
    message (STATUS ${i})  
endforeach()
panoramix_add_library(Panoramix ${panoramix_sources})
target_link_libraries (Panoramix ${DEPENDENCY_LIBS})


# the test project
message(STATUS "panoramix_test_sources:")
foreach(i ${panoramix_test_sources})
    message (STATUS ${i})  
endforeach()
find_package(GTest REQUIRED)
enable_testing ()
include_directories (${DEPENDENCY_INCLUDES})
include_directories (${GTEST_INCLUDE_DIR})
panoramix_add_executable(Panoramix.UnitTest ${panoramix_test_sources} ./panoramix.unittest.hpp ./panoramix.unittest.cpp)
target_link_libraries (Panoramix.UnitTest Panoramix)
target_link_libraries (Panoramix.UnitTest optimized ${GTEST_LIBRARY})
target_link_libraries (Panoramix.UnitTest optimized ${GTEST_MAIN_LIBRARY})
target_link_libraries (Panoramix.UnitTest debug ${GTEST_LIBRARY_DEBUG})
target_link_libraries (Panoramix.UnitTest debug ${GTEST_MAIN_LIBRARY_DEBUG})
target_link_libraries (Panoramix.UnitTest ${DEPENDENCY_LIBS})

