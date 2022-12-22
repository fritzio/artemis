# enable_testing()

# add_executable(hello_test hello_test.cpp)

# target_link_libraries(hello_test GTest::gtest_main)

# include(GoogleTest)
# gtest_discover_tests(hello_test)

set(BINARY ${CMAKE_PROJECT_NAME}_test)

FILE(GLOB_RECURSE TEST_FILES "*.cpp")

add_executable(${BINARY} ${TEST_FILES})

add_test(NAME ${BINARY} COMMAND ${BINARY})

target_link_libraries(${BINARY} GTest::gtest_main)
target_link_libraries(${BINARY} artemis)

include(GoogleTest)
gtest_discover_tests(${BINARY})

# target_link_libraries(${BINARY} PUBLIC gtest)