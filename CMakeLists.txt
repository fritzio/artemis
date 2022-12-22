cmake_minimum_required(VERSION 3.6)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(artemis VERSION 0.0.1 DESCRIPTION "WIP")

# set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -g -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wformat")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wparentheses")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wmultichar")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wtrigraphs")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wpointer-arith")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wreturn-type")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-function")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-narrowing")

include_directories(core)
include_directories(examples)
include_directories(test)

# set the output directory of the binaries to the current folder
# set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

# include_directories(third-party/thrust)
# include_directories(third-party/googletest)

set(CMAKE_PREFIX_PATH third-party/thrust)

find_package(Thrust REQUIRED CONFIG)
# find_package(GoogleTest)

thrust_create_target(Thrust FROM_OPTIONS)

#set(THRUST_DIRECTORY third-party/thrust)
#add_subdirectory(${THRUST_DIRECTORY})

# add_definitions(-D THRUST_HOST_SYSTEM=CPP)
# add_definitions(-D THRUST_DEVICE_SYSTEM=CPP)

# include(cmake/thrust.cmake)

add_subdirectory(core)
add_subdirectory(examples)
add_subdirectory(test)

add_subdirectory(third-party/thrust)
add_subdirectory(third-party/googletest)

# file(GLOB_RECURSE CORE_FILES "core/*.cpp")

# add_library(artemis SHARED ${CORE_FILES})

# set_target_properties(artemis PROPERTIES VERSION ${PROJECT_VERSION})

#target_link_libraries(artemis PUBLIC thrust)

# add_subdirectory(core)
# add_subdirectory(test)
# add_subdirectory(third-party/gtest)
