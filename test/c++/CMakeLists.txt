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

# Using OpenMP
include_directories(/usr/local/include)
find_package(OpenMP REQUIRED)
if(APPLE)
    set(CMAKE_C_COMPILER clang)
    set(CMAKE_CXX_COMPILER clang++)

    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY omp)
    endif()
endif()
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#windows vs2017 doesn't need opm library link 
set(LIBomp "")
if (APPLE)
  set(LIBomp "omp")
endif()

aux_source_directory(. SRC_MAIN)
include_directories(${SWEOS_DIR}/include)
link_directories(${SWEOS_DIR}/lib /usr/local/lib)
add_definitions(-DUSE_OMP=1)
list(APPEND Examples main test_LookUpTableForest test_H2ONaCl_LUT test_omp)
foreach(example IN LISTS Examples)
    message(STATUS "example ${example}")
    add_executable(${example}  ${example}.cpp)
    target_link_libraries(${example} PRIVATE eosH2ONaCl_par ${LIBomp})
endforeach()

message(STATUS "")
message(STATUS "SWEOS_DIR lib path: ${SWEOS_DIR}/lib")
message(STATUS "SWEOS ${SWEOS_VERSION} has been configured for ${SWEOS_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")

