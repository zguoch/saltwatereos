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
set(SWEOS_LIB_PATH "../Library/lib")
set(SWEOS_LIB_NAME "eosH2ONaCl")

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
include_directories(/usr/local/Cellar/libomp/11.0.0/include /usr/local/Cellar/libomp/11.1.0/include)
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
    set_config_option(HAVE_OpenMP "OpenMP")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
#windows vs2017 doesn't need opm library link 
set(LIBomp "")
if (APPLE)
  set(LIBomp "omp")
endif()
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# FreeSteam
set(LIBFreeSteam "")
if (FreeSteam)
    set_config_option(HAVE_FreeSteam "FreeSteam")
    # add freesteam library path
    link_directories("/usr/local/lib")
    set(LIBFreeSteam "freesteam")
endif()
# set icon
set(myApp_ICON ${CMAKE_CURRENT_SOURCE_DIR}/Icons/logo.icns)

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
link_directories(${SWEOS_LIB_PATH} "/usr/local/lib" "${Qt5_DIR}/../../../lib")
FILE(GLOB HDRS_SWEOS_Library "${PROJECT_SOURCE_DIR}/../Library/include/*.h" "${PROJECT_SOURCE_DIR}/../Library/include/*.H")
FILE(GLOB SWEOS_LIB_FILE "${SWEOS_LIB_PATH}/*${SWEOS_LIB_NAME}*")

# Use the include path and library for Qt that is used by VTK.
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set your files and resources here
set( Srcs main.cxx  MainWindow.cxx qroundprogressbar.cpp  Meter.cpp)
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
  WIN32
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
# test: disable macos bundle, cancle it after debug
# set(ENABLE_OS_SPECIFIC_INSTALL ON)
if(APPLE AND NOT ENABLE_OS_SPECIFIC_INSTALL)
  message("Configure for macOS install")
  set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE ON)
  set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)
  set_source_files_properties(${QON_QM_FILES} PROPERTIES MACOSX_PACKAGE_LOCATION "Resources/translations")
  set(MACOSX_BUNDLE_ICON_FILE logo.icns)
  set_source_files_properties(${myApp_ICON} PROPERTIES 
    MACOSX_PACKAGE_LOCATION "Resources")
endif()

target_link_libraries(${PROGRAM_NAME} Qt5::Core Qt5::Gui Qt5::Widgets Qt5::Concurrent ${SWEOS_LIB_NAME} ${LIBFreeSteam} ${LIBomp})
target_link_libraries(${PROGRAM_NAME} ${VTK_LIBRARIES} )


########################################################
###                   INSTALL                        ###
########################################################
if(WIN32)
  if(ENABLE_OS_SPECIFIC_INSTALL)
    set(SWEOS_BIN .)
    set(SWEOS_LIB .)
    set(SWEOS_DOC .)
    set(SWEOS_MAN .)
    set(SWEOS_INC .)
  else()
    include(GNUInstallDirs)
    set(SWEOS_BIN ${CMAKE_INSTALL_BINDIR})
    set(SWEOS_LIB ${CMAKE_INSTALL_LIBDIR})
    set(SWEOS_DOC ${CMAKE_INSTALL_DOCDIR})
    set(SWEOS_MAN ${CMAKE_INSTALL_MANDIR}/man1)
    set(SWEOS_INC ${CMAKE_INSTALL_INCLUDEDIR})
  endif()
  if(CYGWIN)
    unix2dos(SWEOS_API)
    if(ENABLE_PRIVATE_API)
      unix2dos(SWEOS_PRIVATE_API)
    endif()
    unix2dos(WELCOME_FILE)
    unix2dos(SDK_FILE)
    unix2dos(LICENSE_FILE)
    unix2dos(CREDITS_FILE)
    unix2dos(CHANGELOG_FILE)
    unix2dos(TUTORIAL_FILES)
    foreach(DIR ${DEMOS_DIRS})
      file(GLOB DEMO_FILES ${DIR}/?*.*)
      unix2dos(DEMO_FILES)
    endforeach()
  endif()
elseif(APPLE AND ENABLE_OS_SPECIFIC_INSTALL)
  # set these so that the files get installed nicely in the MacOSX
  # .app bundle
  message("macOS install path config")
  set(SWEOS_BIN ../MacOS)
  set(SWEOS_LIB ../MacOS)
  set(SWEOS_DOC ../../..)
  set(SWEOS_MAN ../../..)
  set(SWEOS_INC ../MacOS)
elseif(APPLE AND NOT ENABLE_OS_SPECIFIC_INSTALL)
  # set these so that the files get installed nicely in the MacOSX
  # .app bundle
  set(SWEOS_BIN .)
  set(SWEOS_LIB ./lib)
  set(SWEOS_DOC ./doc)
  set(SWEOS_MAN ./man)
  set(SWEOS_INC ./include)
else()
  include(GNUInstallDirs)
  set(SWEOS_BIN ${CMAKE_INSTALL_BINDIR})
  set(SWEOS_LIB ${CMAKE_INSTALL_LIBDIR})
  set(SWEOS_DOC ${CMAKE_INSTALL_DOCDIR})
  set(SWEOS_MAN ${CMAKE_INSTALL_MANDIR}/man1)
  set(SWEOS_INC ${CMAKE_INSTALL_INCLUDEDIR})
endif() 

install (TARGETS ${PROGRAM_NAME} DESTINATION ${SWEOS_BIN})
# install(FILES ${Hdrs_bash} DESTINATION ${SWEOS_INC})
install(FILES ${HDRS_SWEOS_Library} DESTINATION ${SWEOS_INC})
install(FILES ${SWEOS_LIB_FILE} DESTINATION ${SWEOS_LIB})
install(FILES ${PROJECT_SOURCE_DIR}/../LICENSE DESTINATION .)

# install API
# 0. c++
FILE(GLOB SWEOS_API_cpp "${PROJECT_SOURCE_DIR}/../test/c++/*.*")
install(FILES ${SWEOS_API_cpp} DESTINATION "API/c++/")
# 1. python
FILE(GLOB SWEOS_API_Python "${PROJECT_SOURCE_DIR}/../Library/API/python/pyswEOS/*")
install(FILES ${SWEOS_API_Python} DESTINATION "API/python/")
# 2. tcl
FILE(GLOB SWEOS_API_tcl "${PROJECT_SOURCE_DIR}/../Library/API/tcl/*")
install(FILES ${SWEOS_API_tcl} DESTINATION "API/tcl/")

########################################################
###                   CPACK                          ###
########################################################
set(CPACK_PACKAGE_VENDOR "Zhikui Guo (zguo@geomar.de)")
set(CPACK_PACKAGE_VERSION_MAJOR ${SWEOS_MAJOR_VERSION})
set(CPACK_PACKAGE_VERSION_MINOR ${SWEOS_MINOR_VERSION})
set(CPACK_PACKAGE_VERSION_PATCH ${SWEOS_PATCH_VERSION})
set(CPACK_PACKAGE_DESCRIPTION_FILE ${WELCOME_FILE})
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY
    "swEOS, a desktop application and API for salt-water (H2ONaCl) EOS and thermodynamic properties.")
if(SWEOS_EXTRA_VERSION MATCHES "-git.*") # so that we'll overwrite the archives
  set(CPACK_PACKAGE_FILE_NAME ${PROGRAM_NAME}-git-${SWEOS_OS})
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ${PROGRAM_NAME}-git-source)
else()
  set(CPACK_PACKAGE_FILE_NAME ${PROGRAM_NAME}-${SWEOS_VERSION}-${SWEOS_OS})
  set(CPACK_SOURCE_PACKAGE_FILE_NAME ${PROGRAM_NAME}-${SWEOS_VERSION}-source)
endif()
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${PROGRAM_NAME}")
set(CPACK_RESOURCE_FILE_LICENSE ${LICENSE_FILE})
set(CPACK_RESOURCE_FILE_README ${WELCOME_FILE})
set(CPACK_RESOURCE_FILE_WELCOME ${WELCOME_FILE})
set(CPACK_PACKAGE_EXECUTABLE "${PROGRAM_NAME}")
if(ENABLE_PACKAGE_STRIP)
  set(CPACK_STRIP_FILES TRUE)
else()
  set(CPACK_STRIP_FILES FALSE)
endif()
set(CPACK_SOURCE_IGNORE_FILES "${CMAKE_CURRENT_BINARY_DIR}" )
if (WIN32)
    if (USE_WIX_TOOLSET)
        set(CPACK_GENERATOR "WIX") # this need WiX Tooset installed and a path to candle.exe
    else ()
        set(CPACK_GENERATOR "NSIS") # this needs NSIS installed, and available
    endif ()
    set(CPACK_SOURCE_GENERATOR "ZIP")
elseif ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    set(CPACK_GENERATOR "DragNDrop")
    set(CPACK_BUNDLE_APPLE_CERT_APP "Developer ID Application: Zhikui GUO (AFCJBNSQCL)")
else ()
    set(CPACK_GENERATOR "DEB")
    set(CPACK_SOURCE_GENERATOR "TGZ")
endif ()

include(CPack)

# Configure versions of macro in head file

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
message(STATUS " * Qt5_DIR: " ${Qt5_DIR})
message(STATUS " * VTK version: " ${SWEOS_VTK_VERSION})
message(STATUS " * C compiler: " ${CMAKE_C_COMPILER})
message(STATUS " * C++ compiler: " ${CMAKE_CXX_COMPILER})
message(STATUS " * Install prefix: " ${CMAKE_INSTALL_PREFIX})
message(STATUS "")

# solve font blurry on UI
# set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/Info.plist)