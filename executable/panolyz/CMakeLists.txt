# files
file (GLOB SOURCES "." *.cpp *.hpp)

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
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})
panoramix_add_executable (Panolyz ${SOURCES})
target_link_libraries (Panolyz ${PROJECT_NAME} ${Qt_LIBS})
set_property(TARGET Panolyz PROPERTY FOLDER "${PROJECT_NAME}.Executable")