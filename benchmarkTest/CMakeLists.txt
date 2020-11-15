# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)


set(PROGRAM_NAME swEOS_test)

project(SaltWaterEOS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


# Using OpenMP
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

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(../Library/include /usr/local/Cellar/libomp/10.0.0/include)
link_directories("../Library/lib" "/usr/local/lib" "/Users/zguo/Qt/5.9.9/clang_64/lib")


#2. Source code dir
aux_source_directory(. SRC_MAIN)

# set( Srcs_test ./SWEOSbash.cpp ./MultiProgressBar.cpp)

# set( Hdrs_bash ./SWEOSbash.h ./MultiProgressBar.h)


# add Executable target
add_executable(${PROGRAM_NAME}  main.cpp)

set(LIBomp "")
if (APPLE)
  set(LIBomp "omp")
endif()
target_link_libraries(${PROGRAM_NAME} eosH2ONaCl ${LIBomp} )


