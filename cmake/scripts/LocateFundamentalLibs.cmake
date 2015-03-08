# set Fundamental_INCLUDES, Fundamental_LIBS, Fundamental_PATHS

#find_package (Eigen3 REQUIRED)
find_package (OpenCV REQUIRED)
#find_package (GLPK REQUIRED)
#find_package (GLEW REQUIRED)
# find_package (Armadillo REQUIRED CONFIG)
find_package (MOSEK REQUIRED)
find_package (SuiteSparse REQUIRED)

set (SuiteSparse_USE_LAPACK_BLAS "on")
include(${USE_SuiteSparse})

set (Fundamental_INCLUDES 
	#${EIGEN3_INCLUDE_DIR} 
	${OpenCV_INCLUDE_DIRS} 
	#${GLPK_INCLUDE_DIR} 
	#${GLEW_INCLUDE_DIRS} 
#	${ARMADILLO_INCLUDE_DIRS}
	${MOSEK_INCLUDE_DIR} 
	";${SuiteSparse_INCLUDE_DIRS}"
)

foreach (i ${Fundamental_INCLUDES})
	message (STATUS "Fundamental includes: " ${i})
endforeach ()

set (Fundamental_LIBS 
	${OpenCV_LIBS} 
	#${GLPK_LIBRARY} 
	#${GLEW_LIBRARIES} 
	#${ARMADILLO_LIBRARIES}
	${MOSEK_LIBRARY}
	"${SuiteSparse_LIBRARIES}"
)

foreach (i ${Fundamental_LIBS})
	message (STATUS "Fundamental libs: " ${i})
endforeach ()

# add path to opencv dlls on win
if (DEFINED _OpenCV_LIB_PATH)
	list (APPEND Fundamental_PATH ${_OpenCV_LIB_PATH} ${MOSEK_BIN_DIR})
endif ()

# if (WIN32)
# 	message (WARNING "BLAS or LAPACK dlls may be required by Armadillo in runtime, remember to add them to PATH")
# endif ()

# list (APPEND Fundamental_PATH ${ARMADILLO_LIBRARY_DIRS})

list (APPEND Fundamental_PATH "${SuiteSparse_LIB_DIRS}")