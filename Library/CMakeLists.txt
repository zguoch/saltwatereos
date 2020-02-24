# 0. CMake Minimum version
cmake_minimum_required (VERSION 2.8)

project(H2ONaCl)

#2. Source code dir
aux_source_directory(. SRC_MAIN)


# library
add_subdirectory(src)
# Head files dir
include_directories(src include)

# add freesteam library path
link_directories("/usr/local/lib")

# add Executable target
add_executable(example main.cpp)

# link library
target_link_libraries(example eosH2ONaCl freesteam)


# install path 
install (TARGETS eosH2ONaCl DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)



