# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)

project(test_swEOS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# set compiler flag
add_definitions("-Wno-array-bounds") # turn off this array bound check warning at this moment, this warning accurs at where dim!=3

# 1. set path of SWEOS library, which can be downloaded from https://github.com/zguoch/saltwatereos/releases according to your OS
set(SWEOS_DIR "${PROJECT_SOURCE_DIR}/../../Library" CACHE FILEPATH "Main path of SWEOS library")
if(EXISTS ${SWEOS_DIR})
  message(STATUS "SWEOS head files found: " ${SWEOS_DIR}/include)  
  message(STATUS "SWEOS libraries found: " ${SWEOS_DIR}/lib/libeosH2ONaCl.a )   
else()
  message(STATUS "SWEOS_DIR=${SWEOS_DIR}")
  message(FATAL_ERROR "Please specify path of SWEOS library which contains lib and include paths\n cmake -DSWEOS_DIR=path_of_SWEOS .. ")
endif()

aux_source_directory(. SRC_MAIN)
include_directories(${SWEOS_DIR}/include)
link_directories(${SWEOS_DIR}/lib)

list(APPEND Examples main test_LookUpTableForest test_H2ONaCl_LUT)
foreach(example IN LISTS Examples)
    message(STATUS "example ${example}")
    add_executable(${example}  ${example}.cpp)
    target_link_libraries(${example} PRIVATE eosH2ONaCl)
endforeach()

message(STATUS "SWEOS_DIR lib path: ${SWEOS_DIR}/lib")

