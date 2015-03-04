
# tests
file (GLOB SOURCES "." *.cpp *.hpp)

# google test
find_package (gtest REQUIRED)

include_directories (${GTEST_INCLUDE_DIR})


# define test data location
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

# build test project
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})
add_executable (Tests ${SOURCES})
target_link_libraries (Tests ${PROJECT_NAME} ${Qt_LIBS})
target_link_libraries (Tests optimized ${GTEST_LIBRARY})
target_link_libraries (Tests optimized ${GTEST_MAIN_LIBRARY})
target_link_libraries (Tests debug ${GTEST_LIBRARY_DEBUG})
target_link_libraries (Tests debug ${GTEST_MAIN_LIBRARY_DEBUG})

add_test (NAME Tests COMMAND Tests)
set_tests_properties (Tests PROPERTIES ENVI`RONMENT "${LD_VARNAME}=${LD_PATH}" )
set_property(TARGET Tests PROPERTY FOLDER "${PROJECT_NAME}.Executable")