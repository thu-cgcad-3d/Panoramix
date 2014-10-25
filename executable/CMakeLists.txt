

# tests
file (GLOB TESTCPPS "." test_*.cpp)


# apps
file (GLOB APPS "." app_*)


# google test
set (GTEST_DIR "gtest-1.7.0")
include_directories (${GTEST_DIR}/include)

# define test data location
get_filename_component(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/" REALPATH)
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

# build test projects
foreach (TCPP ${TESTCPPS})
	get_filename_component (T ${TCPP} NAME_WE)
	message (STATUS ${T})
	include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})
	add_executable (${T} ${T}.cpp)
	target_link_libraries (${T} ${CORE_STATIC_LIB} ${DERIV_STATIC_LIB} ${VIS_STATIC_LIB} ${REC_STATIC_LIB} ${Qt_LIBS})
	add_test (NAME ${T} COMMAND ${T})
	set_tests_properties (${T} PROPERTIES ENVIRONMENT "${LD_VARNAME}=${LD_PATH}" )
	set_property(TARGET ${T} PROPERTY FOLDER "${PROJECT_NAME}.Tests")
endforeach()

# build app projects
foreach (TFILE ${APPS})
	get_filename_component (T ${TFILE} NAME_WE)
	message (STATUS ${T})
	include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})
	if (IS_DIRECTORY "${TFILE}")
		file (GLOB TSOURCES 
			"${TFILE}/*.cpp" 
			"${TFILE}/*.hpp" 
			"${TFILE}/*.cc" 
			"${TFILE}/*.c" 
			"${TFILE}/*.hh" 
			"${TFILE}/*.h")
		message (STATUS "sources: ${TSOURCES}")
		add_executable (${T} ${TSOURCES})
	else ()
		add_executable (${T} ${T}.cpp)
	endif ()
	target_link_libraries (${T} ${CORE_STATIC_LIB} ${DERIV_STATIC_LIB} ${VIS_STATIC_LIB} ${REC_STATIC_LIB} ${Qt_LIBS})
	set_property(TARGET ${T} PROPERTY FOLDER "${PROJECT_NAME}.Apps")
endforeach()

