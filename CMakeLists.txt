# 0. CMake Minimum version
cmake_minimum_required (VERSION 2.8)

# 1. Project Information
project (Demo)
# Version of project
set (Demo_VERSION_MAJOR 1)
set (Demo_VERSION_MINOR 0)

# 2. Executable target
#add_executable(Demo0 main.cpp)
#add_executable(Demo1 main1.cpp func1.cpp)

#2. Source code dir
aux_source_directory(. SRC_MAIN)


# library
add_subdirectory(library)
# Head files dir
include_directories(library)

# add freesteam library path
link_directories("/usr/local/lib")

# add Executable target
add_executable(example main.cpp)

# link library
target_link_libraries(example eosH2ONaCl freesteam)




