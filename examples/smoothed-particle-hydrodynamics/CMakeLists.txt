FILE(GLOB_RECURSE SPH_FILES "*.cpp")

add_executable(smoothed-particle-hydrodynamics ${SPH_FILES})

target_link_libraries(smoothed-particle-hydrodynamics PUBLIC artemis)