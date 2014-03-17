set(src_dir "@CMAKE_CURRENT_SOURCE_DIR@")
configure_file (
	"${src_dir}/config.hpp.cmake" 
	"${src_dir}/config.hpp"
)

file (GLOB_RECURSE HEADERS ./*.hpp)
file (GLOB_RECURSE SOURCES ./*.cpp)

foreach (H ${HEADERS})
	message (STATUS "Header: " ${H})
endforeach ()

foreach (S ${SOURCES})
	message (STATUS "Source: " ${S})
endforeach ()


add_library (${PROJECT_LIB_STATIC_NAME} STATIC ${SOURCES} ${HEADERS})
#add_library (${PROJECT_LIB_SHARED_NAME} SHARED ${SOURCES})

target_link_libraries (${PROJECT_LIB_STATIC_NAME} ${THIRDPARTY_LIBS})
#target_link_libraries (${PROJECT_LIB_SHARED_NAME} ${THIRDPARTY_LIBS})

