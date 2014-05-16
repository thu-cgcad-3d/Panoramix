set (GTEST_DIR "gtest-1.7.0")
include_directories (${GTEST_DIR}/include)

# set all tests
set (TESTNAMES 
	test_basic
	test_util
	test_optimization
	test_template_utilities
	test_expression
	test_mesh 
	test_feature
	test_views_net 
	test_views_net_incremented
	test_views_net_calibrate_cameras
	test_views_net_outdoor
	test_regions_net
	test_vis2d
	test_vis3d
)


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
foreach (T ${TESTNAMES})
	include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})
	add_executable (${T} ${T}.cpp)
	target_link_libraries (${T} ${CORE_STATIC_LIB} ${DERIV_STATIC_LIB} ${VIS_STATIC_LIB} ${REC_STATIC_LIB} ${Qt_LIBS})
	add_test (NAME ${T} COMMAND ${T})
	set_tests_properties (${T} PROPERTIES ENVIRONMENT "${LD_VARNAME}=${LD_PATH}" )
	set_property(TARGET ${T} PROPERTY FOLDER "${PROJECT_NAME}.Tests")
endforeach()


# sandbox TEST
set (SBTESTNAME test_sandbox)
include_directories (${Fundamental_INCLUDES} ${Qt_INCLUDES})
add_executable(${SBTESTNAME} ${SBTESTNAME}.cpp)
target_link_libraries (${SBTESTNAME} ${SANDBOX_STATIC_LIB} ${Qt_LIBS})
add_test (NAME ${SBTESTNAME} COMMAND ${SBTESTNAME})
set_tests_properties (${SBTESTNAME} PROPERTIES ENVIRONMENT "${LD_VARNAME}=${LD_PATH}" )
set_property(TARGET ${SBTESTNAME} PROPERTY FOLDER "${PROJECT_NAME}.Tests")


