# set Fundamental_INCLUDES, Fundamental_LIBS, Fundamental_PATHS

find_package (Eigen3 REQUIRED)
find_package (OpenCV REQUIRED)
find_package (GLPK REQUIRED)

set (Fundamental_INCLUDES ${EIGEN3_INCLUDE_DIR} ${OpenCV_INCLUDE_DIRS} ${GLPK_INCLUDE_DIR} )
foreach (i ${Fundamental_INCLUDES})
	message (STATUS "Fundamental includes: " ${i})
endforeach ()
set (Fundamental_LIBS ${OpenCV_LIBS} ${GLPK_LIBRARY})
foreach (i ${Fundamental_LIBS})
	message (STATUS "Fundamental libs: " ${i})
endforeach ()
# add path to opencv dlls on win
if (DEFINED _OpenCV_LIB_PATH)
	set (Fundamental_PATH ${_OpenCV_LIB_PATH})
endif ()