set(THRUST_HOST_SYSTEM "CPP" CACHE STRING "Selects the host system")
set_property(CACHE THRUST_HOST_SYSTEM PROPERTY STRINGS "CPP" "TBB" "OMP")
add_definitions(-D THRUST_HOST_SYSTEM=${THRUST_HOST_SYSTEM})
message(STATUS "THRUST_HOST_SYSTEM : ${THRUST_HOST_SYSTEM}")

set(THRUST_DEVICE_SYSTEM "CPP" CACHE STRING "Selects the device system")
set_property(CACHE THRUST_DEVICE_SYSTEM PROPERTY STRINGS "CUDA" "CPP" "TBB" "OMP")
add_definitions(-D THRUST_DEVICE_SYSTEM=${THRUST_DEVICE_SYSTEM})
message(STATUS "THRUST_DEVICE_SYSTEM : ${THRUST_DEVICE_SYSTEM}")
