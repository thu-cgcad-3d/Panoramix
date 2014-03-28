set (GTEST_DIR "gtest-1.7.0")
include_directories (${GTEST_DIR}/include)

# set all tests
set (TESTNAMES test_basic)
set (TESTNAMES_WITH_VIS test_util test_mesh test_feature test_views_net test_regions_net)

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
	message (WARNING "CMake cannot configure Visual Studio to modify the environment path during program execution, "
		"therefore you have to do this manually: add 'PATH=\$(PATH);${LD_PATH};' to [Project Property]->[Debug]->[Environment] ")
endif ()


foreach (T ${TESTNAMES})
	include_directories (${Fundamental_INCLUDES})
	add_executable (${T} ${T}.cpp)
	target_link_libraries (${T} ${CORE_STATIC_LIB})
	add_test (NAME ${T} COMMAND ${T})
	set_tests_properties (${T} PROPERTIES ENVIRONMENT "${LD_VARNAME}=${LD_PATH}" )
	# input vars and environment vars do not apply on IDEs (VS and XCode), manual settings are needed
endforeach ()


foreach (T ${TESTNAMES_WITH_VIS})
	include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})
	add_executable (${T} ${T}.cpp)
	target_link_libraries (${T} ${CORE_STATIC_LIB} ${VIS_STATIC_LIB} ${Qt_LIBS})
	add_test (NAME ${T} COMMAND ${T})
	set_tests_properties (${T} PROPERTIES ENVIRONMENT "${LD_VARNAME}=${LD_PATH}" )
endforeach()

