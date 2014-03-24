
set (SRC_MODULES core)

foreach (M ${SRC_MODULES})
	message (STATUS "module name: " ${M})
	file (GLOB_RECURSE M_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.hpp)
	file (GLOB_RECURSE M_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${M}/*.cpp)
	foreach (H ${M_HEADERS})
		message (STATUS "	header: " ${H})
	endforeach ()
	foreach (S ${M_SOURCES})
		message (STATUS "	source: " ${S})
	endforeach ()
	list (APPEND HEADERS ${M_HEADERS})
	list (APPEND SOURCES ${M_SOURCES})
	source_group (${M} FILES ${M_HEADERS} ${M_SOURCES})
endforeach ()


add_library (${PROJECT_LIB_STATIC_NAME} STATIC ${SOURCES} ${HEADERS})
#add_library (${PROJECT_LIB_SHARED_NAME} SHARED ${SOURCES})

target_link_libraries (${PROJECT_LIB_STATIC_NAME} ${THIRDPARTY_LIBS})
#target_link_libraries (${PROJECT_LIB_SHARED_NAME} ${THIRDPARTY_LIBS})
