cmake_minimum_required(VERSION 3.10)

# Setup the project.
# Its name is defined to be the controller directory name.
set(PROJECT "lqr")
project(${PROJECT})

set( CMAKE_DISABLE_IN_SOURCE_BUILD ON)
set( CMAKE_CXX_STANDARD 14)
set( CMAKE_C_STANDARD 11)
set( EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}) 
set( LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/build/lib) 

# # Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message("Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for test_footStest_footPlacementtepPlan-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

# Get C or C++ sources in the current directory (only).
file(GLOB C_SOURCES ./src/*.c)
file(GLOB CPP_SOURCES ./src/*.cpp)
set(SOURCES ${C_SOURCES} ${CPP_SOURCES})

# Link with the Webots controller library.
link_directories($ENV{WEBOTS_HOME}/lib/controller)
set (LIBRARIES ${CMAKE_SHARED_LIBRARY_PREFIX}Controller${CMAKE_SHARED_LIBRARY_SUFFIX} ${CMAKE_SHARED_LIBRARY_PREFIX}CppController${CMAKE_SHARED_LIBRARY_SUFFIX})
include_directories(
$ENV{WEBOTS_HOME}/include/controller/c 
$ENV{WEBOTS_HOME}/include/controller/cpp)

# Setup the target executable.
add_executable(${PROJECT} ${SOURCES})
target_link_libraries(${PROJECT} ${LIBRARIES})

# add_custom_command(TARGET ${PROJECT} POST_BUILD 
# COMMAND ${CMAKE_COMMAND} -E
#   copy ${PROJECT_SOURCE_DIR}/build/bin/* ${PROJECT_SOURCE_DIR})