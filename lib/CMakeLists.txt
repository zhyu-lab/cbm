cmake_minimum_required(VERSION 3.0)

# Build the matrix library
include_directories(matrix)

# Build the split library
include_directories(split)
add_library(split split/split.cpp)

# Build the config library
include_directories(config)
add_library(config config/Config.cpp)

# Build the threadpool library
include_directories(threadpool)
add_library(threadpool threadpool/ThreadPool.cpp)

# Build the mydefine library
include_directories(mydefine)
add_library(mydefine mydefine/MyDefine.cpp)
target_link_libraries(mydefine config clonecaller inputparser threadpool)

# Build the inputparser library
include_directories(inputparser)
add_library(inputparser inputparser/InputParser.cpp)
target_link_libraries(inputparser mydefine)

# Build the clonecaller library
include_directories(clonecaller)
add_library(clonecaller clonecaller/CloneCaller.cpp)
target_link_libraries(clonecaller mydefine split)
