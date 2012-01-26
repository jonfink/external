cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
set(ROS_BUILD_TYPE Release)

rosbuild_init()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()

IF(APPLE)
  SET(SOCP_LIBS ${SOCP_LIBS} "-framework Accelerate")  # or "-framework accelerate" ?
  MESSAGE(STATUS "MacOS X detected. Added '-framework Accelerate' to compiler flags")
ELSE(APPLE)
  SET(SOCP_LIBS ${SOCP_LIBS} "-lcblas")
ENDIF(APPLE)

#common commands for building c++ executables and libraries
rosbuild_add_library(${PROJECT_NAME} src/socp.c)
target_link_libraries(${PROJECT_NAME} ${SOCP_LIBS})
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
rosbuild_add_executable(example src/example.c)
target_link_libraries(example ${PROJECT_NAME})

rosbuild_add_library(armadillo_socp src/armadillo_socp.cc src/socp.c)
target_link_libraries(armadillo_socp ${SOCP_LIBS})

rosbuild_add_executable(armadillo_example src/armadillo_example.cc)
target_link_libraries(armadillo_example armadillo_socp)

rosbuild_add_executable(radio_example src/radio_example.cc)
target_link_libraries(radio_example armadillo_socp)

rosbuild_add_executable(radio_example2 src/radio_example2.cc)
target_link_libraries(radio_example2 armadillo_socp)
