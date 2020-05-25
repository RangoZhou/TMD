set (BUILD_MODE NORMAL)
cmake_minimum_required (VERSION 3.5)

project(TMD)
# The version number.
set (TMD_VERSION_MAJOR 1)
set (TMD_VERSION_MINOR 0)
# configure a header file to pass some of the CMake settings to the source code
configure_file (
  ${PROJECT_SOURCE_DIR}/src/TMDConfig.cmake.template
  ${PROJECT_SOURCE_DIR}/src/TMDConfig.h
)


set(BOOST_ROOT "${PROJECT_SOURCE_DIR}/lib/boost_1_71_0_installation")
set(BOOST_INCLUDE "${BOOST_ROOT}/include")
set(BOOST_LIB "${BOOST_ROOT}/lib")
set(LINKED_LIBS "${BOOST_LIB}/libboost_filesystem.a" "${BOOST_LIB}/libboost_system.a" "${BOOST_LIB}/libboost_thread.a" "${BOOST_LIB}/libboost_serialization.a" "${BOOST_LIB}/libboost_program_options.a" "-lstdc++fs" "-lpthread")
# message("${LINKED_LIBS}")

set(TMD_SOURCE_INCLUDE "${PROJECT_SOURCE_DIR}/src")
set(TMD_SOURCE_DIR "${PROJECT_SOURCE_DIR}/src")

set(NR_INCLUDE "${PROJECT_SOURCE_DIR}/src/nr")

# add the binary tree to the search path for include files
include_directories("${NR_INCLUDE}")
include_directories("${TMD_SOURCE_INCLUDE}")
include_directories("${BOOST_INCLUDE}")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# collect all source files
# AUX_SOURCE_DIRECTORY(${TMD_SOURCE_DIR} SOURCE_FILES)

set(TMD_SRCS ${TMD_SOURCE_DIR}/vec3d.cpp)
# set(TMD_SRCS ${TMD_SOURCE_DIR}/vec3d.cpp ${TMD_SOURCE_DIR}/score.cpp)

set(TMD_SOURCE_FILES ${TMD_SOURCE_DIR}/main.cpp ${TMD_SRCS})

# add the executable
# link library
if(${BUILD_MODE} STREQUAL "NORMAL")
set(CMAKE_CXX_FLAGS "-O3 -std=c++17")
add_executable(tmd ${TMD_SOURCE_FILES})
target_link_libraries(tmd ${LINKED_LIBS})
endif()

if(${BUILD_MODE} STREQUAL "DEBUG")
set(CMAKE_CXX_FLAGS "-g -std=c++17")
add_executable(tmd_debug ${TMD_SOURCE_FILES})
target_link_libraries(tmd_debug ${LINKED_LIBS})
endif()

if(${BUILD_MODE} STREQUAL "PROFILE")
set(CMAKE_CXX_FLAGS "-O3 -std=c++17")
add_executable(tmd_profile ${TMD_SOURCE_FILES})
target_link_libraries(tmd_profile ${LINKED_LIBS} "-lprofiler")
endif()