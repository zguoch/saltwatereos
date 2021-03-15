# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)

project(test_swEOS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. set path of SWEOS library, which can be downloaded from https://github.com/zguoch/saltwatereos/releases according to your OS
set(SWEOS_DIR "../../" CACHE FILEPATH "Main path of SWEOS library")
if(EXISTS ${SWEOS_DIR})
  message(STATUS "SWEOS head files found: " ${SWEOS_DIR}/include)  
  message(STATUS "SWEOS libraries found: " ${SWEOS_DIR}/lib/libeosH2ONaCl.a )   
else()
  message(FATAL_ERROR "Please specify path of SWEOS library which contains lib and include paths\n cmake -DSWEOS_DIR=path_of_SWEOS .. ")
endif()


aux_source_directory(. SRC_MAIN)
include_directories(${SWEOS_DIR}/include)
link_directories(${SWEOS_DIR}/lib)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} eosH2ONaCl)
