
add_subdirectory ("gpc")
add_subdirectory ("android-camera")

set (ThirdParty_INCLUDES 
	${CMAKE_CURRENT_SOURCE_DIR}/cereal/include
	${CMAKE_CURRENT_SOURCE_DIR}/gpc
	${CMAKE_CURRENT_SOURCE_DIR}/android-camera/feature_mos/src
	${CMAKE_CURRENT_SOURCE_DIR}/android-camera/feature_stab/db_vlvm
	${CMAKE_CURRENT_SOURCE_DIR}/dlib
	${CMAKE_CURRENT_SOURCE_DIR}/misc PARENT_SCOPE)

set (ThirdParty_LIBS "gpc" "android-camera" PARENT_SCOPE)