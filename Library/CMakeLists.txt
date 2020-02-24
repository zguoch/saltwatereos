# 0. CMake Minimum version
cmake_minimum_required (VERSION 2.8)

project(H2ONaCl)

#2. Source code dir
aux_source_directory(. SRC_MAIN)


# library
add_subdirectory(src)
# Head files dir
include_directories(src include)



# add Executable target
add_executable(example main.cpp)

option(USE_FreeSteam "whether use freesteam library" OFF)
if (USE_FreeSteam)
    # add freesteam library path
    link_directories("/usr/local/lib")
    # example target link library
    target_link_libraries(example eosH2ONaCl freesteam)
else (USE_FreeSteam)
    # example target link library
    target_link_libraries(example eosH2ONaCl)
endif()



# install path 
install (TARGETS eosH2ONaCl DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/lib)



