file(GLOB_RECURSE CORE_FILES "*.cpp")

add_library(artemis SHARED ${CORE_FILES})

target_link_libraries(artemis Thrust)

set_target_properties(artemis PROPERTIES VERSION ${PROJECT_VERSION})