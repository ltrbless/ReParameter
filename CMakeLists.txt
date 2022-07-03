cmake_minimum_required(VERSION 3.0)

aux_source_directory(./src sources)

add_executable(main ${sources} main.cpp)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/libigl/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/cli11)

target_link_libraries(main -lpthread)

