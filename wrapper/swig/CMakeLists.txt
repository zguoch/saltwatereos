# 0. CMake Minimum version
cmake_minimum_required (VERSION 3.12)
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

# 1. SWIG 
FIND_PACKAGE(SWIG REQUIRED)  
INCLUDE(${SWIG_USE_FILE})  
# 1.1 find python
FIND_PACKAGE(PythonLibs)  
SET(PYTHON_INCLUDE_PATH "/Users/zguo/.pyenv/versions/anaconda3-2019.10/include/python3.7m")
SET(PYTHON_LIBRARIES "/Users/zguo/.pyenv/versions/anaconda3-2019.10/lib/libpython3.7m.dylib")
if(PYTHONLIBS_FOUND)
    message("Python include found" ${PYTHON_INCLUDE_PATH} )
    message("Python libs found" ${PYTHON_LIBRARIES} )
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
    INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
    SET(CMAKE_SWIG_FLAGS "")
endif(PYTHONLIBS_FOUND)

#2. Source code dir
aux_source_directory(. SRC_MAIN)
include_directories(${PROJECT_SOURCE_DIR})

# SET_SOURCE_FILES_PROPERTIES(example.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(example.i PROPERTIES SWIG_FLAGS "-includeall")
set ( CMAKE_C_FLAGS "-flat_namespace") # important
# SWIG_ADD_LIBRARY(example python example.i example.c)
swig_add_library(example
        LANGUAGE python
        SOURCES example.i example.c)

SWIG_LINK_LIBRARIES(example ${PYTHON_LIBRARIES})

# install
# execute_process(COMMAND python -c "from distutils.sysconfig import get_python_lib; print get_python_lib()" OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE)
message("python site: " ${PYTHON_SITE_PACKAGES} ${CMAKE_BINARY_DIR})
install(TARGETS example DESTINATION ${PROJECT_SOURCE_DIR}/lib)
install(FILES ${CMAKE_BINARY_DIR}/example.py DESTINATION ${PROJECT_SOURCE_DIR}/lib)

# # install path 
# FILE(GLOB SC_FILES "${PROJECT_SOURCE_DIR}/example.c")
# add_library (example ${SC_FILES})
# # add_library(Fluid STATIC ${PROJECT_SOURCE_DIR}/../../Library/src/Fluid.C)
# # target_link_libraries(example Fluid )
# install (TARGETS example DESTINATION ${PROJECT_SOURCE_DIR}/lib)



