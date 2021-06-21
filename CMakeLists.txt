cmake_minimum_required(VERSION 3.6)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(artemis VERSION 0.0.1 DESCRIPTION "WIP")

set(CMAKE_CXX_STANDARD 17)

include_directories(core)
include_directories(third-party)

include_directories(third-party/thrust)

set(CMAKE_PREFIX_PATH third-party/thrust)

find_package(Thrust REQUIRED CONFIG)

# thrust_create_target(Thrust HOST CPP DEVICE OMP)
thrust_create_target(Thrust FROM_OPTIONS)


#set(THRUST_DIRECTORY third-party/thrust)
#add_subdirectory(${THRUST_DIRECTORY})

# add_definitions(-D THRUST_HOST_SYSTEM=CPP)
# add_definitions(-D THRUST_DEVICE_SYSTEM=CPP)

# include(cmake/thrust.cmake)

add_subdirectory(core)
add_subdirectory(examples)

# file(GLOB_RECURSE CORE_FILES "core/*.cpp")

# add_library(artemis SHARED ${CORE_FILES})

# set_target_properties(artemis PROPERTIES VERSION ${PROJECT_VERSION})

#target_link_libraries(artemis PUBLIC thrust)

# add_subdirectory(core)
# add_subdirectory(test)
# add_subdirectory(third-party/gtest)
