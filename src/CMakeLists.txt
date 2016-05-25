set_property(GLOBAL PROPERTY USE_FOLDERS On)

set (modules core gui misc experimental)
set(CMAKE_AUTOMOC on)
set(CMAKE_AUTOUIC on)
set(CMAKE_AUTORCC on)

set (panoramix_sources "")
set (panoramix_test_sources "")
foreach (M ${modules})
    file (GLOB sources 
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.cpp" 
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.cu" 
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.cuh"
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.qrc"
    )
    file (GLOB test_sources 
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.test.cpp" 
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.test.hpp"
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.test.cu" 
        "${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.test.cuh"
    )
    if (test_sources)
        list (REMOVE_ITEM sources ${test_sources})
    endif ()
    source_group ("Modules\\${M}" FILES ${sources})
    source_group ("Modules\\${M}" FILES ${test_sources})
    list (APPEND panoramix_sources ${sources})
    list (APPEND panoramix_test_sources ${test_sources})
endforeach ()


# the lib project
message(STATUS "panoramix_sources:")
foreach(i ${panoramix_sources})
    message (STATUS ${i})  
endforeach()
panoramix_add_library(Panoramix ${panoramix_sources})
target_include_directories(Panoramix PUBLIC ${DEPENDENCY_INCLUDES})
target_link_libraries (Panoramix ${DEPENDENCY_LIBS})
add_dependencies(Panoramix ${DEPENDENCY_NAMES})


# the test project
message(STATUS "panoramix_test_sources:")
foreach(i ${panoramix_test_sources})
    message (STATUS ${i})  
endforeach()
enable_testing ()
panoramix_add_executable(Panoramix.UnitTest ${panoramix_test_sources} 
    ./panoramix.unittest.hpp ./panoramix.unittest.cpp)
target_include_directories(Panoramix.UnitTest PUBLIC ${DEPENDENCY_INCLUDES})
target_link_libraries (Panoramix.UnitTest Panoramix)
target_link_libraries (Panoramix.UnitTest ${DEPENDENCY_LIBS})
add_dependencies(Panoramix.UnitTest Panoramix)
add_dependencies(Panoramix.UnitTest ${DEPENDENCY_NAMES})

