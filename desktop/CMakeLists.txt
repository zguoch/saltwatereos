cmake_minimum_required(VERSION 3.3...3.12 FATAL_ERROR)
foreach(p
    CMP0071 # 3.10: Let AUTOMOC and AUTOUIC process GENERATED files
    )
  if(POLICY ${p})
    cmake_policy(SET ${p} NEW)
  endif()
endforeach()

project(swEOS LANGUAGES CXX)
set(PROGRAM_NAME ${PROJECT_NAME})
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)


if(${CMAKE_VERSION} VERSION_LESS "3.3") 
else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()

macro(opt OPTION HELP VALUE)
  option(ENABLE_${OPTION} ${HELP} ${VALUE})
  set(OPT_TEXI "${OPT_TEXI}\n@item ENABLE_${OPTION}\n${HELP} (default: ${VALUE})")
endmacro()
macro(set_config_option VARNAME STRING)
  set(${VARNAME} TRUE)
  list(APPEND CONFIG_OPTIONS ${STRING})
  message(STATUS "Found " ${STRING})
endmacro()
opt(OpenMP "Enable OpenMP to parallel computing" ON)
option(FreeSteam "whether use freesteam library" OFF)

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

if(NOT DATE)
  string(TIMESTAMP DATE "%Y-%m-%d")
endif()
set(SWEOS_DATE "${DATE}")

if(NOT SWEOS_HOST)
  execute_process(COMMAND hostname OUTPUT_VARIABLE HOSTNAME
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT HOSTNAME)
    set(HOSTNAME "unknown")
  endif()
  set(SWEOS_HOST "${HOSTNAME}")
endif()

if(NOT SWEOS_PACKAGER)
  execute_process(COMMAND whoami OUTPUT_VARIABLE PACKAGER
                  OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT PACKAGER)
    set(PACKAGER "unknown")
  endif()
  string(REPLACE "\\" " " PACKAGER ${PACKAGER})
  set(SWEOS_PACKAGER "${PACKAGER}")
endif()

if(APPLE)
  set(SWEOS_OS "MacOSX")
elseif(CYGWIN OR MSYS)
  # detect if we use the MinGW compilers on Cygwin - if we do, handle the build
  # as a pure Windows build
  if(CMAKE_CXX_COMPILER_ID MATCHES "GNU" OR
     CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    execute_process(COMMAND ${CMAKE_CXX_COMPILER} -dumpmachine
                    OUTPUT_VARIABLE CXX_COMPILER_MACHINE
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CXX_COMPILER_MACHINE MATCHES "mingw")
      set(SWEOS_OS "Windows")
      set(WIN32 1)
      add_definitions(-DWIN32)
    endif()
  endif()
else()
  set(SWEOS_OS "${CMAKE_SYSTEM_NAME}")
endif()

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
    set_config_option(HAVE_OpenMP "OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
include_directories(/usr/local/Cellar/libomp/11.0.0/include /usr/local/Cellar/libomp/11.1.0/include)

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
message("VTK version: " "${VTK_VERSION}")
message("VTK QT version: " "${VTK_QT_VERSION}")
include(${VTK_USE_FILE})
if("${VTK_VERSION}" STREQUAL "")
  message(FATAL_ERROR "Could not find VTK")
else()
  set(SWEOS_VTK_VERSION ${VTK_VERSION})
endif()
if("${VTK_QT_VERSION}" STREQUAL "")
  message(FATAL_ERROR "VTK was not built with Qt")
else()
  set(SWEOS_QT_VERSION ${VTK_QT_VERSION})
endif()
# add H2ONaCl source and library 
include_directories(../Library/include ../commandline)
link_directories("../Library/lib" "/usr/local/lib" "${Qt5_DIR}/../../../lib")



# Use the include path and library for Qt that is used by VTK.
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set your files and resources here
# qroundprogressbar.cpp
set( Srcs main.cxx  MainWindow.cxx qroundprogressbar.cpp  Meter.cpp)
# qroundprogressbar.h
set( Hdrs MainWindow.h qroundprogressbar.h  Meter.h)

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

set(LIBFreeSteam "")
 
if (FreeSteam)
    set_config_option(HAVE_FreeSteam "FreeSteam")
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/swEOSConfig.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/swEOSConfig.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../swEOSVersion.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/swEOSVersion.h)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../swEOSVersion.h.in
               ${CMAKE_CURRENT_SOURCE_DIR}/../commandline/swEOSVersion.h)

message(STATUS "")
message(STATUS "SWEOS ${SWEOS_VERSION} has been configured for ${SWEOS_OS}")
message(STATUS "")
message(STATUS " * Build type: " ${CMAKE_BUILD_TYPE})
message(STATUS " * Qt version: " ${SWEOS_QT_VERSION})
message(STATUS " * VTK version: " ${SWEOS_VTK_VERSION})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")

# solve font blurry on UI
set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)