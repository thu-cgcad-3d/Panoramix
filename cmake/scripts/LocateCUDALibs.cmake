

if (${USE_CUDA})
find_package(CUDA REQUIRED)
endif ()

macro(panoramix_add_executable PROJECT_NAME)
	if (${USE_CUDA})
		include_directories (${CUDA_TOOLKIT_INCLUDE})
		cuda_add_executable(${PROJECT_NAME} ${ARGN})
	else()
		add_executable(${PROJECT_NAME} ${ARGN})
	endif()
endmacro()

macro(panoramix_add_library PROJECT_NAME)
	if (${USE_CUDA})
		include_directories (${CUDA_TOOLKIT_INCLUDE})
		cuda_add_library(${PROJECT_NAME} ${ARGN})
	else()
		add_library(${PROJECT_NAME} ${ARGN})
	endif()
endmacro()
