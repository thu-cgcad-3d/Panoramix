

# set all tests
set (TESTNAMES test_config test_mesh test_feature)

# define test data location
get_filename_component(TEST_DATA_DIR "${CMAKE_CURRENT_SOURCE_DIR}/data/" REALPATH)
message (STATUS "test data directory: " ${TEST_DATA_DIR})
if (MSVC)
	add_definitions ( "/DPROJECT_TEST_DATA_DIR_STR=\"${TEST_DATA_DIR}\"" )
else ()
	add_definitions ( "-DPROJECT_TEST_DATA_DIR_STR=\"${TEST_DATA_DIR}\"" )
endif ()


# the JOIN FUNCTION
# JOIN("${somelist}" ":" output)
function(JOIN VALUES GLUE OUTPUT)
  string (REGEX REPLACE "([^\\]|^);" "\\1${GLUE}" _TMP_STR "${VALUES}")
  #string (REGEX REPLACE "[\\](.)" "\\1" _TMP_STR "${_TMP_STR}") #fixes escaping
  set (${OUTPUT} "${_TMP_STR}" PARENT_SCOPE)
endfunction()

# set env PATH for dlls/dylibs etc ...
if (WIN32)
	message (STATUS "Win32")
	set (LD_VARNAME "PATH")
	JOIN("${THIRDPARTY_PATH}" ";" LD_PATH)
	set (LD_PATH "$ENV{PATH};${LD_PATH}")
elseif (APPLE)
	message (STATUS "Apple")
	set (LD_VARNAME "DYLD_LIBRARY_PATH")
	JOIN("${THIRDPARTY_PATH}" ":" LD_PATH)
	set (LD_PATH "$ENV{DYLD_LIBRARY_PATH}:${LD_PATH}")
elseif (UNIX)
	message (STATUS "Unix other than Apple")
	set (LD_VARNAME "LD_LIBRARY_PATH")
	JOIN("${THIRDPARTY_PATH}" ":" LD_PATH)
	set (LD_PATH "$ENV{LD_LIBRARY_PATH}:${LD_PATH}")
endif ()
message (STATUS "test ${LD_VARNAME} = ${LD_PATH}")


foreach (T ${TESTNAMES})
	add_executable (${T} ${T}.cpp)
	target_link_libraries (${T} ${PROJECT_LIB_STATIC_NAME})
	add_test (NAME ${T} COMMAND ${T})
	set_tests_properties (${T} PROPERTIES ENVIRONMENT "${LD_VARNAME}=${LD_PATH}" )
	# input vars and environment vars do not apply on IDEs (VS and XCode), manual settings are needed
endforeach ()

