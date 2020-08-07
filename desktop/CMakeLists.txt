cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)
foreach(p
    CMP0071 # 3.10: Let AUTOMOC and AUTOUIC process GENERATED files
    )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

set(PROGRAM_NAME swEOS)

project(SaltWaterEOS LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Using OpenMP
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

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
include_directories(/usr/local/Cellar/libomp/10.0.0/include)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

# set icon
set(MACOSX_BUNDLE_ICON_FILE logo.icns)
set(myApp_ICON ${CMAKE_CURRENT_SOURCE_DIR}/Icons/logo.icns)
set_source_files_properties(${myApp_ICON} PROPERTIES 
    MACOSX_PACKAGE_LOCATION "Resources")




find_package(VTK COMPONENTS
  vtkCommonCore
  vtkFiltersCore
  vtkInfovisCore
  vtkInteractionStyle
  vtkViewsQt
  vtkChartsCore
  vtkCommonDataModel
  vtkRenderingContext2D
  vtkRenderingContextOpenGL2
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingGL2PSOpenGL2
  vtkRenderingOpenGL2
  vtkViewsContext2D
  vtkIOLegacy
)
include(${VTK_USE_FILE})

# add H2ONaCl source and library 
include_directories(../Library/include ../commandline)
link_directories("../Library/lib" "/usr/local/lib" "/Users/zguo/Qt/5.9.9/clang_64/lib")

if("${VTK_QT_VERSION}" STREQUAL "")
  message(FATAL_ERROR "VTK was not built with Qt")
endif()

# Use the include path and library for Qt that is used by VTK.
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set your files and resources here

set( Srcs main.cxx  MainWindow.cxx qroundprogressbar.cpp Meter.cpp)

set( Hdrs MainWindow.h qroundprogressbar.h Meter.h)

set( MOC_Hdrs ui_MainWindow.h )

set( UIs MainWindow.ui )

set( QRCs Icons/icons.qrc)

set( Srcs_bash ../commandline/SWEOSbash.cpp ../commandline/MultiProgressBar.cpp)

set( Hdrs_bash ../commandline/SWEOSbash.h ../commandline/MultiProgressBar.h)


# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# The rest should just work (sure...)
# We have ui files, this will bring in the macro: qt5_wrap_ui
find_package(Qt5 COMPONENTS Widgets Core Concurrent LinguistTools REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})
qt5_add_resources(QRC_Srcs ${QRCs} )

#languages
set(TS_FILES languages/zh_CN.ts)
qt5_add_translation(QON_QM_FILES ${TS_FILES})
set_source_files_properties(${QON_QM_FILES} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources/languages")

source_group("Resources" FILES
  ${UIs}
  ${QRCs}
  ${EXE_ICON} # Not present
)

source_group("Generated" FILES
  ${UI_Srcs}
  ${MOC_Srcs}
  ${QRC_Srcs}
)


# if want to generage mac os app, add MACOSX_BUNDLE follow ${PROGRAM_NAME} 
add_executable(
  ${PROGRAM_NAME}  
  MACOSX_BUNDLE
  ${Srcs} 
  ${Hdrs} 
  ${Srcs_bash} 
  ${Hdrs_bash} 
  ${UI_Srcs} 
  ${MOC_Hdrs} 
  ${QRC_Srcs}
  ${myApp_ICON}
  ${QON_QM_FILES}
  )
option(USE_FreeSteam "whether use freesteam library" OFF)
set(LIBFreeSteam "")
 
if (USE_FreeSteam)
    # add freesteam library path
    link_directories("/usr/local/lib")
    set(LIBFreeSteam "freesteam")
endif()

#windows vs2017 doesn't need opm library link 
set(LIBomp "")
if (APPLE)
  set(LIBomp "omp")
endif()
target_link_libraries(${PROGRAM_NAME} Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Concurrent eosH2ONaCl ${LIBFreeSteam} ${LIBomp})
target_link_libraries(${PROGRAM_NAME} ${VTK_LIBRARIES} )



# solve font blurry on UI
set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)