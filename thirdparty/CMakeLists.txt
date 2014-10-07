
add_subdirectory ("gpc")
add_subdirectory ("android-camera")
add_subdirectory ("gco-v3.0")
add_subdirectory ("lsd_1.6")
add_subdirectory ("slic")

set (ThirdParty_INCLUDES 
	${CMAKE_CURRENT_SOURCE_DIR}/eigen
	${CMAKE_CURRENT_SOURCE_DIR}/cereal/include
	${CMAKE_CURRENT_SOURCE_DIR}/gpc
	${CMAKE_CURRENT_SOURCE_DIR}/gco-v3.0
	${CMAKE_CURRENT_SOURCE_DIR}/lsd_1.6
	${CMAKE_CURRENT_SOURCE_DIR}/slic
	${CMAKE_CURRENT_SOURCE_DIR}/android-camera/feature_mos/src
	${CMAKE_CURRENT_SOURCE_DIR}/android-camera/feature_stab/db_vlvm
	${CMAKE_CURRENT_SOURCE_DIR}/dlib
	${CMAKE_CURRENT_SOURCE_DIR}/misc PARENT_SCOPE)

set (ThirdParty_LIBS "gpc" "android-camera" "gco-v3.0" "lsd_1.6" "slic" PARENT_SCOPE)
