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

set(CMAKE_INCLUDE_CURRENT_DIR ON)

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
include_directories(../src)
link_directories("../lib" "/usr/local/lib")

if("${VTK_QT_VERSION}" STREQUAL "")
  message(FATAL_ERROR "VTK was not built with Qt")
endif()

# Use the include path and library for Qt that is used by VTK.
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# Set your files and resources here
set( Srcs main.cxx MainWindow.cxx )

set( Hdrs MainWindow.h )

set( MOC_Hdrs ui_MainWindow.h )

set( UIs MainWindow.ui )

set( QRCs Icons/icons.qrc )


# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# The rest should just work (sure...)
# We have ui files, this will bring in the macro: qt5_wrap_ui
find_package(Qt5 COMPONENTS Widgets REQUIRED QUIET)
qt5_wrap_ui(UI_Srcs ${UIs})
qt5_add_resources(QRC_Srcs ${QRCs} )

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
# 如果不设置MACOSX_BUNDLE， 窗口变化后VTK渲染会卡顿. 但是设置了之后字体会变粗，估计是与系统设置有关
add_executable(${PROGRAM_NAME}  MACOSX_BUNDLE
  ${Srcs} ${Hdrs} ${UI_Srcs} ${MOC_Hdrs} ${QRC_Srcs})
target_link_libraries(${PROGRAM_NAME} Qt5::Core Qt5::Gui Qt5::Widgets eosH2ONaCl freesteam)
target_link_libraries(${PROGRAM_NAME} ${VTK_LIBRARIES})

# solve font blurry on UI
set_target_properties(${PROGRAM_NAME} PROPERTIES MACOSX_BUNDLE_INFO_PLIST "/Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/Info.plist")