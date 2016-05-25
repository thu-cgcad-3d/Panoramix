add_subdirectory ("lsd_1.6")
add_subdirectory ("slic")
add_subdirectory ("vanishingpoint-msac")
add_subdirectory ("vpdetection")

set (ThirdParty_INCLUDES 
	${CMAKE_CURRENT_SOURCE_DIR}/eigen
	${CMAKE_CURRENT_SOURCE_DIR}/cereal/include
	${CMAKE_CURRENT_SOURCE_DIR}/lsd_1.6
	${CMAKE_CURRENT_SOURCE_DIR}/slic
	${CMAKE_CURRENT_SOURCE_DIR}/vanishingpoint-msac
	${CMAKE_CURRENT_SOURCE_DIR}/vpdetection
	${CMAKE_CURRENT_SOURCE_DIR}/misc 
	PARENT_SCOPE)

set (ThirdParty_LIBS 
	"lsd_1.6" 
	"slic" 
	"vanishingpoint-msac" 
	"vpdetection" 
	PARENT_SCOPE
)
