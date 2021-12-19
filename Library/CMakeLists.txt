# 0. CMake Minimum version
cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)
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

project(eosH2ONaCl LANGUAGES CXX C)
# project(${PROJECT_NAME} LANGUAGES CXX)
set(LIB_NAME "${PROJECT_NAME}")
if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

# options
option(USE_OMP "Enable OpenMP to parallel computing" OFF)

# Using OpenMP
if(USE_OMP)
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
  add_definitions(-DUSE_OMP=1)
  set(LIB_NAME "${PROJECT_NAME}_par")
else()
  set(LIB_NAME "${PROJECT_NAME}")
endif()

# set compiler flag
if (APPLE)
  add_definitions("-Wno-array-bounds -Wno-unused-command-line-argument") # turn off this array bound check warning at this moment, this warning accurs at where dim!=3
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
# add_subdirectory(swig/H2O)
# add_subdirectory(swig/NaCl)
# add_subdirectory(swig/H2ONaCl) # split omp-related code to another file, and do not compile that file in swig

# install path 
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
link_directories(/usr/local/lib)
FILE(GLOB SC_FILES "${PROJECT_SOURCE_DIR}/src/PROST/*.c" "${PROJECT_SOURCE_DIR}/src/PROST/*.h" "${PROJECT_SOURCE_DIR}/src/*.H" "${PROJECT_SOURCE_DIR}/src/*.cpp")

# default is STATIC library, for Mex compilation on Linux system, the SHARED library is required!!! But static lib works for MAC OS
# compile both STATIC and SHARED library
list(APPEND buildTypes STATIC SHARED)
foreach(buildType IN LISTS buildTypes)
  add_library ("${LIB_NAME}_${buildType}" ${buildType} ${SC_FILES} )
  set_target_properties("${LIB_NAME}_${buildType}" PROPERTIES OUTPUT_NAME ${LIB_NAME})
  target_link_libraries("${LIB_NAME}_${buildType}" ${LIBomp})
  install (TARGETS "${LIB_NAME}_${buildType}" DESTINATION ${PROJECT_SOURCE_DIR}/lib)
endforeach()

# app
list(APPEND Apps lutInfo lutLookup lut2vtu lutGen_2D lutGen_3D lutRandom)
foreach(app IN LISTS Apps)
    message(STATUS "Add app ${app}")
    add_executable(${app}  ${SC_FILES} "app/${app}.cpp")
    target_link_libraries(${app} ${LIBomp})
    install (TARGETS ${app} DESTINATION ${PROJECT_SOURCE_DIR}/bin)
endforeach()
install (FILES "${PROJECT_SOURCE_DIR}/app/randSample.py" DESTINATION ${PROJECT_SOURCE_DIR}/bin)

# testing
enable_testing()
add_executable(test_lut ${SC_FILES} test/test_H2ONaCl_LUT.cpp)
target_link_libraries(test_lut ${LIBomp})
add_test(test_lut0 test_lut 1)

message(STATUS "")
message(STATUS "SWEOS ${SWEOS_VERSION} has been configured for ${SWEOS_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
if(USE_OMP)
  message(STATUS " * Use OpenMP: " ${LIBomp})
endif(USE_OMP)
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")
# ===================================================================================================


