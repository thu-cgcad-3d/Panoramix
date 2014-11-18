set (MATLAB_INCLUDES "")
set (MATLAB_LIBS "")

if (${USE_MATLAB})
find_package(MATLAB REQUIRED)
#  MATLAB_INCLUDE_DIR: include path for mex.h, engine.h
#  MATLAB_LIBRARIES:   required libraries: libmex, etc
#  MATLAB_MEX_LIBRARY: path to libmex.lib
#  MATLAB_MX_LIBRARY:  path to libmx.lib
#  MATLAB_MAT_LIBRARY:  path to libmat.lib # added
#  MATLAB_ENG_LIBRARY: path to libeng.lib
#  MATLAB_ROOT: path to Matlab's root directory
list (APPEND MATLAB_INCLUDES ${MATLAB_INCLUDE_DIR})
list (APPEND MATLAB_LIBS ${MATLAB_LIBRARIES})
endif ()

foreach (i ${MATLAB_INCLUDES})
	message (STATUS "MATLAB includes: " ${i})
endforeach ()
foreach (i ${MATLAB_LIBS})
	message (STATUS "MATLAB libs: " ${i})
endforeach ()
