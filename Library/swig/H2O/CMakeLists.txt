
# ================== whether build wrapper for other language, e.g. python, tcl =====================
option(Build_API_MultiLanguage "whether build API for other langulages" ON)
if(Build_API_MultiLanguage)
    FIND_PACKAGE(SWIG REQUIRED)  
    INCLUDE(${SWIG_USE_FILE})  
    FIND_PACKAGE(PythonLibs)  
    if(PYTHONLIBS_FOUND)
        message("Python include found" ${PYTHON_INCLUDE_PATH} )
        message("Python libs found" ${PYTHON_LIBRARIES} )
        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
        SET(CMAKE_SWIG_FLAGS "")
    endif(PYTHONLIBS_FOUND)
    include_directories("${PROJECT_SOURCE_DIR}/include")
    SET_SOURCE_FILES_PROPERTIES(H2O.i PROPERTIES CPLUSPLUS ON)
    SET_SOURCE_FILES_PROPERTIES(H2O.i PROPERTIES SWIG_FLAGS "-includeall")
    # set ( CMAKE_C_FLAGS "-flat_namespace") # important
    set (CMAKE_CXX_FLAGS "-flat_namespace") # important
    # ------ 1. python ------------------
    option(Build_API_Python "whether build API for Python" ON)
    if(Build_API_Python)
        swig_add_library(H2O
        LANGUAGE python
        SOURCES H2O.i ${PROJECT_SOURCE_DIR}/src/H2O.C ${PROJECT_SOURCE_DIR}/src/Fluid.C )
        SWIG_LINK_LIBRARIES(H2O ${PYTHON_LIBRARIES})
        # install
        install(TARGETS H2O DESTINATION ${PROJECT_SOURCE_DIR}/API/python)
        install(FILES ${CMAKE_CURRENT_BINARY_DIR}/H2O.py DESTINATION ${PROJECT_SOURCE_DIR}/API/python)
        install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/test.py DESTINATION ${PROJECT_SOURCE_DIR}/API/python)
    endif(Build_API_Python)
    # ------ 2. tcl ------------------
    option(Build_API_tcl "whether build API for tcl" ON)
    if(Build_API_tcl)
        swig_add_library(tclH2O
        LANGUAGE tcl
        SOURCES H2O.i ${PROJECT_SOURCE_DIR}/src/H2O.C ${PROJECT_SOURCE_DIR}/src/Fluid.C )
        SWIG_LINK_LIBRARIES(tclH2O tcl)
        # install
        set_target_properties( tclH2O PROPERTIES OUTPUT_NAME "H2O" )
        install(TARGETS tclH2O DESTINATION ${PROJECT_SOURCE_DIR}/API/tcl RENAME "H2O")
        install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/test.tcl DESTINATION ${PROJECT_SOURCE_DIR}/API/tcl)
    endif(Build_API_tcl)
    # ------------------------------------
endif(Build_API_MultiLanguage)

# ===================================================================================================


