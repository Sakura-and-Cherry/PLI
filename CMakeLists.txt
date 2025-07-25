cmake_minimum_required(VERSION 3.16)
project(PGM)


# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mbmi2")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w") # Disable all warnings
set(CMAKE_CXX_STANDARD 17)

# sdsl-lite
include_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/include)
link_directories(${PROJECT_SOURCE_DIR}/include/sdsl-lite/lib)

include_directories(${PROJECT_SOURCE_DIR}/include)


# pgm_build
add_executable(pgm_build_uint32 pgm_build.cpp)
target_compile_definitions(pgm_build_uint32 PRIVATE USE_UINT32_KEY)
target_link_libraries(pgm_build_uint32 sdsl)

add_executable(pgm_build_uint64 pgm_build.cpp)
target_compile_definitions(pgm_build_uint64 PRIVATE USE_UINT64_KEY)
target_link_libraries(pgm_build_uint64 sdsl)

# data_prepare
add_executable(data_prepare_uint32 data_prepare.cpp)
target_compile_definitions(data_prepare_uint32 PRIVATE USE_UINT32_KEY)

add_executable(data_prepare_uint64 data_prepare.cpp)
target_compile_definitions(data_prepare_uint64 PRIVATE USE_UINT64_KEY)

# pgm_query
add_executable(pgm_query_uint32 pgm_query.cpp)
target_compile_definitions(pgm_query_uint32 PRIVATE USE_UINT32_KEY)
target_link_libraries(pgm_query_uint32 sdsl)

add_executable(pgm_query_uint64 pgm_query.cpp)
target_compile_definitions(pgm_query_uint64 PRIVATE USE_UINT64_KEY)
target_link_libraries(pgm_query_uint64 sdsl)


add_executable(test test.cpp)
target_link_libraries(test sdsl)