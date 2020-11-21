# 0. CMake Minimum version
cmake_minimum_required (VERSION 2.12)
if (POLICY CMP0086)
  # cmake warns if loaded from a min-3.0-required parent dir, so silence the warning:
  cmake_policy(SET CMP0086 NEW)
endif()
if (POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()

project(H2ONaCl)
project(SaltWaterEOS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#2. Source code dir
aux_source_directory(. SRC_MAIN)

# source
add_subdirectory(src)
add_subdirectory(swig/H2O)
add_subdirectory(swig/NaCl)
add_subdirectory(swig/H2ONaCl)
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


# target_include_directories(water PUBLIC include)

# install path 
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
FILE(GLOB SC_FILES "${PROJECT_SOURCE_DIR}/src/PROST/*.c" "${PROJECT_SOURCE_DIR}/src/PROST/*.h" "${PROJECT_SOURCE_DIR}/src/*.H" "${PROJECT_SOURCE_DIR}/src/*.cpp")
add_library (eosH2ONaCl ${SC_FILES} src/H2O.cpp src/Fluid.cpp src/NaCl.cpp)
# add_library(water STATIC src/H2O.cpp)
# target_link_libraries(eosH2ONaCl water Fluid NaCl)
install (TARGETS eosH2ONaCl DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# # ===================================================================================================


