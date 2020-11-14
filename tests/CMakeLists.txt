cmake_minimum_required(VERSION 3.10)

if (NOT TARGET gtest)
  include(googletest)
  __fetch_googletest(
    ${PROJECT_SOURCE_DIR}/cmake
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}_googletest
)
endif()

add_executable(
  ${PROJECT_NAME}_unit_tests
  TestHungarianAlgorithm.cpp
)

target_link_libraries(
  ${PROJECT_NAME}_unit_tests
  PUBLIC
    ${LIBRARY_NAME}
  PRIVATE
    gtest_main
)

add_test(
  NAME
    ${PROJECT_NAME}_unit_tests
  COMMAND
    $<TARGET_FILE:${PROJECT_NAME}_unit_tests>
)
