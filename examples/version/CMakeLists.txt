FILE(GLOB_RECURSE VERSION_FILES "*.cpp")
set(BINARY version)

add_executable(${BINARY} ${VERSION_FILES})

target_link_libraries(${BINARY} PUBLIC artemis)