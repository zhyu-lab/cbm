cmake_minimum_required(VERSION 3.0)
project(CBM)

SET(CMAKE_BUILD_TYPE Release)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -std=c++11")

add_subdirectory(lib)

include_directories (
	${CBM_SOURCE_DIR}/lib/config
	${CBM_SOURCE_DIR}/lib/clonecaller
	${CBM_SOURCE_DIR}/lib/inputparser
	${CBM_SOURCE_DIR}/lib/matrix
	${CBM_SOURCE_DIR}/lib/mydefine
	${CBM_SOURCE_DIR}/lib/split
	${CBM_SOURCE_DIR}/lib/threadpool
)

add_executable(cbm src/cbm.cpp)

target_link_libraries(cbm mydefine)

SET_TARGET_PROPERTIES(cbm PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CBM_BINARY_DIR}/bin)

find_package(Threads REQUIRED)
if(THREADS_HAVE_PTHREAD_ARG)
	set_property(TARGET cbm PROPERTY COMPILE_OPTIONS "-pthread")
	set_property(TARGET cbm PROPERTY INTERFACE_COMPILE_OPTIONS "-pthread")
endif()
if(CMAKE_THREAD_LIBS_INIT)
	target_link_libraries(cbm "${CMAKE_THREAD_LIBS_INIT}")
endif()

INSTALL(TARGETS cbm
	RUNTIME DESTINATION bin
)
