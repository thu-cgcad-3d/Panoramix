# files

if (MSVC)
	add_definitions ( "/DPROJECT_TEST_DATA_DIR_STR=\"${TEST_DATA_DIR}\"" )
	add_definitions(/DMATLAB_MEX_FILE) #define matlab macros
	add_definitions(/DMX_COMPAT_32)
else ()
	add_definitions ( "-DPROJECT_TEST_DATA_DIR_STR=\"${TEST_DATA_DIR}\"" )
	add_definitions(-DMATLAB_MEX_FILE) #define matlab macros
	add_definitions(-DMX_COMPAT_32)
endif ()

JOIN("${Fundamental_PATH}" ";" LD_PATH)
  configure_file (${CMAKE_CURRENT_SOURCE_DIR}/setpath.m.in ${MATLAB_CODE_DIR}/setpath.m)

macro(panoramix_add_mex MEXNAME MFILE)
  include_directories (${Fundamental_INCLUDES} ${MATLAB_INCLUDES})
  add_library(${MEXNAME} SHARED ../Matlabdef.def ${ARGN})
  target_link_libraries(${MEXNAME} 
    ${PROJECT_NAME} 
    ${Fundamental_LIBS} 
    ${MATLAB_LIBS})
  set_property(GLOBAL PROPERTY USE_FOLDERS On)
  set_property(TARGET ${MEXNAME} PROPERTY FOLDER "${PROJECT_NAME}.Mex")

  # 32-bit or 64-bit mex
  if(WIN32)
    if (CMAKE_CL_64)
        SET_TARGET_PROPERTIES(${MEXNAME} PROPERTIES SUFFIX .mexw64)
    else(CMAKE_CL_64)
        SET_TARGET_PROPERTIES(${MEXNAME} PROPERTIES SUFFIX .mexw32)
    endif(CMAKE_CL_64)
  else(WIN32)
    if (CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${MEXNAME} PROPERTIES SUFFIX .mexa64 PREFIX "")
    else(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${MEXNAME} PROPERTIES SUFFIX .mexglx PREFIX "")
    endif (CMAKE_SIZEOF_VOID_P MATCHES "8")
  endif(WIN32)

  # install to /bin by default
  install(TARGETS ${MEXNAME} DESTINATION ${MATLAB_CODE_DIR})
  install(FILES ${MFILE} DESTINATION ${MATLAB_CODE_DIR})
endmacro()

add_subdirectory (boundary_graph)
add_subdirectory (mixed_graph)