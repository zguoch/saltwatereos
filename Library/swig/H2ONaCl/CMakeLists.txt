
# ================== whether build wrapper for other language, e.g. python, tcl =====================
option(Build_API_MultiLanguage "whether build API for other langulages" ON)
set(specie "H2ONaCl")
INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/src/PROST)
file(GLOB PROST_SRC "${PROJECT_SOURCE_DIR}/src/PROST/*.c")
file(GLOB SWEOS_SRC "${PROJECT_SOURCE_DIR}/src/*.cpp")
set(API_SRC ${PROST_SRC} ${SWEOS_SRC})
if(Build_API_MultiLanguage)
    FIND_PACKAGE(SWIG REQUIRED)  
    INCLUDE(${SWIG_USE_FILE})  
    include_directories("${PROJECT_SOURCE_DIR}/include")
    SET_SOURCE_FILES_PROPERTIES(${specie}.i PROPERTIES CPLUSPLUS ON)
    SET_SOURCE_FILES_PROPERTIES(${specie}.i PROPERTIES SWIG_FLAGS "-includeall")
    # set ( CMAKE_C_FLAGS "-flat_namespace") # important
    set (CMAKE_CXX_FLAGS "-flat_namespace") # important
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
        swig_add_library(${specie}
        LANGUAGE python
        SOURCES ${specie}.i ${API_SRC})
        SWIG_LINK_LIBRARIES(${specie} ${PYTHON_LIBRARIES})
        # install
        install(TARGETS ${specie} DESTINATION ${PROJECT_SOURCE_DIR}/API/python)
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${specie}.py DESTINATION ${PROJECT_SOURCE_DIR}/API/python)
        install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/test_${specie}.py DESTINATION ${PROJECT_SOURCE_DIR}/API/python)
    endif(Build_API_Python)
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

