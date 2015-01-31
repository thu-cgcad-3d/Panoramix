# files
file (GLOB SOURCES "." *.cpp *.hpp)

# define test data location
get_filename_component(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../tests/data/" REALPATH)
message (STATUS "test data directory: " ${TEST_DATA_DIR})
if (MSVC)
	add_definitions ( "/DPROJECT_TEST_DATA_DIR_STR=\"${TEST_DATA_DIR}\"" )
else ()
	add_definitions ( "-DPROJECT_TEST_DATA_DIR_STR=\"${TEST_DATA_DIR}\"" )
endif ()

if (MSVC)
	JOIN("${Fundamental_PATH};${Qt_PATH}" ";" LD_PATH)
	configure_file (${CMAKE_CURRENT_SOURCE_DIR}/setpath.bat.in ${CMAKE_CURRENT_BINARY_DIR}/setpath.bat)
	message (WARNING "CMake cannot configure Visual Studio to modify the environment path during program execution, "
		"therefore you have to do this manually: add 'PATH=\$(PATH);${LD_PATH};' to [Project Property]->[Debug]->[Environment] ")
endif ()


set_property(GLOBAL PROPERTY USE_FOLDERS On)

# build project
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES} "../../src" "../../tests")
add_executable (SceneAnalyzer ${SOURCES})
target_link_libraries (SceneAnalyzer ${PROJECT_NAME} ${Qt_LIBS})
set_property(TARGET SceneAnalyzer PROPERTY FOLDER "${PROJECT_NAME}.Executable")