# sandbox module

include_directories (${Fundamental_INCLUDES})

file (GLOB SANDBOX_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
file (GLOB SANDBOX_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)

source_group ("sandbox" FILES ${SANDBOX_HEADERS} ${SANDBOX_SOURCES})

add_library (${SANDBOX_STATIC_LIB} STATIC ${SANDBOX_SOURCES} ${SANDBOX_HEADERS})
target_link_libraries (${SANDBOX_STATIC_LIB} ${Fundamental_LIBS})

