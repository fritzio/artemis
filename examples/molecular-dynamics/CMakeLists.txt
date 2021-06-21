FILE(GLOB_RECURSE MD_FILES "*.cpp")

add_executable(molecular-dynamics ${MD_FILES})

target_link_libraries(molecular-dynamics PUBLIC artemis)