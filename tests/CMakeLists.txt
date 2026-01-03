add_executable(unit_tests core_test.cpp)
# GTest::gtest_main をリンクすると、main関数を書かずに済みます
target_link_libraries(unit_tests PRIVATE GTest::gtest_main mphhc_lib)

include(GoogleTest)
gtest_discover_tests(unit_tests)
