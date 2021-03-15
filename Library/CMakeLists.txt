# 0. CMake Minimum version
cmake_minimum_required (VERSION 2.12)
if (POLICY CMP0086)
  # cmake warns if loaded from a min-3.0-required parent dir, so silence the warning:
  cmake_policy(SET CMP0086 NEW)
endif()
if (POLICY CMP0078)
  cmake_policy(SET CMP0078 NEW)
endif()
if (POLICY CMP0042)
  cmake_policy(SET CMP0042 NEW)
endif()

project(eosH2ONaCl CXX C)
# project(${PROJECT_NAME} LANGUAGES CXX)
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

set(SWEOS_MAJOR_VERSION 1)
set(SWEOS_MINOR_VERSION 7)
set(SWEOS_PATCH_VERSION 0)
set(SWEOS_EXTRA_VERSION "")
set(SWEOS_EXTRA_VERSION_TEXI "${SWEOS_EXTRA_VERSION}")
if(NOT SWEOS_RELEASE)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} log -1 --format=%h
                    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR} ERROR_QUIET
                    OUTPUT_VARIABLE GIT_COMMIT_HASH
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  if(GIT_COMMIT_HASH)
    set(SWEOS_EXTRA_VERSION "${SWEOS_EXTRA_VERSION}-git-${GIT_COMMIT_HASH}")
  endif()
  set(SWEOS_EXTRA_VERSION_TEXI "${SWEOS_EXTRA_VERSION_TEXI} (development version)")
endif()
set(SWEOS_VERSION "${SWEOS_MAJOR_VERSION}.${SWEOS_MINOR_VERSION}")
set(SWEOS_VERSION "${SWEOS_VERSION}.${SWEOS_PATCH_VERSION}${SWEOS_EXTRA_VERSION}")
set(SWEOS_SHORT_LICENSE "GNU General Public License")

#2. Source code dir
aux_source_directory(. SRC_MAIN)
# Head files dir
include_directories(src include)
# source
add_subdirectory(src)
add_subdirectory(swig/H2O)
add_subdirectory(swig/NaCl)
add_subdirectory(swig/H2ONaCl)

# install path 
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
FILE(GLOB SC_FILES "${PROJECT_SOURCE_DIR}/src/PROST/*.c" "${PROJECT_SOURCE_DIR}/src/PROST/*.h" "${PROJECT_SOURCE_DIR}/src/*.H" "${PROJECT_SOURCE_DIR}/src/*.cpp")

add_library (${PROJECT_NAME}  ${SC_FILES} )
# add_library(water STATIC src/H2O.cpp)
target_link_libraries(${PROJECT_NAME})
install (TARGETS ${PROJECT_NAME} DESTINATION ${PROJECT_SOURCE_DIR}/lib)

message(STATUS "")
message(STATUS "SWEOS ${SWEOS_VERSION} has been configured for ${SWEOS_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
# ===================================================================================================


