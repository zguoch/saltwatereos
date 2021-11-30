# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)


set(PROGRAM_NAME examples)

project(examples LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# p4est lib
set(P4EST_ROOT "../../ThirdParties/p4est-2.8/local")
set(P4EST_INCLUDE_DIR "${P4EST_ROOT}/include")
set(P4EST_LIB_DIR "${P4EST_ROOT}/lib")
list(APPEND P4EST_LIBS p4est sc)

# SWEOS lib
set(SWEOS_DIR "../../Library")
set(SWEOS_INCLUDE_DIR "${SWEOS_DIR}/include")
set(SWEOS_LIB_DIR "${SWEOS_DIR}/lib")
list(APPEND SWEOS_LIBS eosH2ONaCl)

#2. add p4est include and lib dir
include_directories(${P4EST_INCLUDE_DIR})
link_directories(${P4EST_LIB_DIR})
# add sweos include and lib
include_directories(${SWEOS_INCLUDE_DIR})
link_directories(${SWEOS_LIB_DIR})
aux_source_directory(. SRC_MAIN)

list(APPEND Examples step1 step2 step3 step4 step5 step6)
foreach(example IN LISTS Examples)
    message(STATUS "example ${example}")
    add_executable(${example}  ${example}.cpp)
    target_link_libraries(${example} ${P4EST_LIBS} ${SWEOS_LIBS})
endforeach()

message(STATUS "")
message(STATUS ${PROGRAM_NAME} " has been configured for ${SWEOS_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
