
add_subdirectory ("gpc")

set (ThirdParty_INCLUDES ${CMAKE_CURRENT_SOURCE_DIR}/cereal/include
	${CMAKE_CURRENT_SOURCE_DIR}/gpc
	${CMAKE_CURRENT_SOURCE_DIR}/dlib
	${CMAKE_CURRENT_SOURCE_DIR}/misc PARENT_SCOPE)

set (ThirdParty_LIBS "gpc" PARENT_SCOPE)