cmake_minimum_required(VERSION 3.10)

add_executable(hungarian_algorithm_app
  ${CMAKE_CURRENT_LIST_DIR}/App.cpp
)

target_link_libraries(hungarian_algorithm_app
  PUBLIC
    ${LIBRARY_NAME}
)
