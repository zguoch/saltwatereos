# 0. CMake Minimum version
cmake_minimum_required (VERSION 2.8)

project(H2ONaCl)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
FILE(GLOB SC_FILES "${PROJECT_SOURCE_DIR}/src/PROST/*.c" "${PROJECT_SOURCE_DIR}/src/PROST/*.h" "${PROJECT_SOURCE_DIR}/src/*.H" "${PROJECT_SOURCE_DIR}/src/*.C")
add_library (eosH2ONaCl ${SC_FILES})
install (TARGETS eosH2ONaCl DESTINATION ${PROJECT_SOURCE_DIR}/lib)



