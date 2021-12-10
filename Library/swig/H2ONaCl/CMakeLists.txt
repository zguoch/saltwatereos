# Using OpenMP
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

# ================== whether build wrapper for other language, e.g. python, tcl =====================
option(Build_API_MultiLanguage "whether build API for other langulages" ON)
set(specie "H2ONaCl")
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/PROST)
file(GLOB PROST_SRC "${PROJECT_SOURCE_DIR}/src/PROST/*.c")
file(GLOB SWEOS_SRC "${PROJECT_SOURCE_DIR}/src/H2ONaCl.cpp" "${PROJECT_SOURCE_DIR}/src/H2O.cpp" "${PROJECT_SOURCE_DIR}/src/NaCl.cpp" "${PROJECT_SOURCE_DIR}/src/Polynomial.cpp" "${PROJECT_SOURCE_DIR}/src/PolynomialRootFinder.cpp")
set(API_SRC ${PROST_SRC} ${SWEOS_SRC})
if(Build_API_MultiLanguage)
    FIND_PACKAGE(SWIG REQUIRED) 
    INCLUDE(${SWIG_USE_FILE})  
    include_directories("${PROJECT_SOURCE_DIR}/include")
    SET_SOURCE_FILES_PROPERTIES(${specie}.i PROPERTIES CPLUSPLUS ON)
    SET_SOURCE_FILES_PROPERTIES(${specie}.i PROPERTIES SWIG_FLAGS "-includeall")
    if (APPLE)
        # set ( CMAKE_C_FLAGS "-flat_namespace") # important
        set (CMAKE_CXX_FLAGS "-flat_namespace") # important
    endif()
    # ------ 1. python ------------------
    option(Build_API_Python "whether build API for Python" ON)
    if(Build_API_Python)
        FIND_PACKAGE(PythonLibs)  
        if(PYTHONLIBS_FOUND)
            message("Python include found: " ${PYTHON_INCLUDE_PATH} )
            message("Python libs found: " ${PYTHON_LIBRARIES} )
            INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
            SET(CMAKE_SWIG_FLAGS "")
        endif(PYTHONLIBS_FOUND)
        link_directories(/usr/local/lib)
        swig_add_library(${specie}
        LANGUAGE python
        SOURCES ${specie}.i ${API_SRC})
        SWIG_LINK_LIBRARIES(${specie} ${PYTHON_LIBRARIES} ${LIBomp})

        # install 
        set(INSTALL_PYTHON_API_PYTHON ${PROJECT_SOURCE_DIR}/API/python/pyswEOS)
        install(TARGETS ${specie} DESTINATION ${INSTALL_PYTHON_API_PYTHON})
        file(GLOB files_install
            "${CMAKE_CURRENT_SOURCE_DIR}/*.py"
            "${CMAKE_CURRENT_SOURCE_DIR}/*.pvsm"
        )
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${specie}.py DESTINATION ${INSTALL_PYTHON_API_PYTHON})
        install(FILES ${files_install} DESTINATION ${INSTALL_PYTHON_API_PYTHON})
    endif(Build_API_Python)
    # --------- Matlab ----------------
    option(Build_API_Matlab "whether build API for Matlab" ON)
    if(Build_API_Matlab)
        # install 
        set(INSTALL_MATLAB_API ${PROJECT_SOURCE_DIR}/API/Matlab)
        file(GLOB files_install_matlab
            "${PROJECT_SOURCE_DIR}/Matlab/${specie}/*.cpp"
            "${PROJECT_SOURCE_DIR}/Matlab/${specie}/*.m"
        )
        install(FILES ${files_install_matlab} DESTINATION ${INSTALL_MATLAB_API}/${specie})
    endif(Build_API_Matlab)
    # ------ 2. tcl ------------------
    option(Build_API_tcl "whether build API for tcl" ON)
    if(Build_API_tcl)
        FIND_PACKAGE(Tcl)  
        if(TCL_FOUND)
            message("TCL include found: " ${TCL_INCLUDE_PATH} )
            message("TCL libs found: " ${TCL_LIBRARY} )
            INCLUDE_DIRECTORIES(${TCL_INCLUDE_PATH})
            SET(CMAKE_SWIG_FLAGS "")
        endif(TCL_FOUND)
        swig_add_library(tcl${specie}
        LANGUAGE tcl
        SOURCES ${specie}.i ${API_SRC})
        SWIG_LINK_LIBRARIES(tcl${specie} ${TCL_LIBRARY})
        # install
        set_target_properties( tcl${specie} PROPERTIES OUTPUT_NAME "${specie}" SUFFIX "")
        install(TARGETS tcl${specie} DESTINATION ${PROJECT_SOURCE_DIR}/API/tcl RENAME "${specie}")
        install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/test_${specie}.tcl DESTINATION ${PROJECT_SOURCE_DIR}/API/tcl)
    endif(Build_API_tcl)
    # ------------------------------------
endif(Build_API_MultiLanguage)

# ===================================================================================================


