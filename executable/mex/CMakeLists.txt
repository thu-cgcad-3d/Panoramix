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

set (MEXNames 
  "mexGetSegmentationGraph"
  "mexGetBoundaryLabels" 
  "mexGetBoundaryFeas")

foreach (CPP_FILE ${MEXNames})
  
  include_directories (${Fundamental_INCLUDES} ${MATLAB_INCLUDES})
  add_library(${CPP_FILE} SHARED ${CPP_FILE}.cpp ./Matlabdef.def)
  target_link_libraries(${CPP_FILE} 
  	${PROJECT_NAME}	
  	${Fundamental_LIBS} 
  	${MATLAB_LIBS})
  set_property(GLOBAL PROPERTY USE_FOLDERS On)
  set_property(TARGET ${CPP_FILE} PROPERTY FOLDER "${PROJECT_NAME}.Executable.Mex")

  # 32-bit or 64-bit mex
  if(WIN32)
    if (CMAKE_CL_64)
        SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexw64)
    else(CMAKE_CL_64)
        SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexw32)
    endif(CMAKE_CL_64)
  else(WIN32)
    if (CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexa64 PREFIX "")
    else(CMAKE_SIZEOF_VOID_P MATCHES "8")
        SET_TARGET_PROPERTIES(${CPP_FILE} PROPERTIES SUFFIX .mexglx PREFIX "")
    endif (CMAKE_SIZEOF_VOID_P MATCHES "8")
  endif(WIN32)

  # install to /bin by default
  install(TARGETS ${CPP_FILE} DESTINATION ${MATLAB_CODE_DIR})

endforeach()