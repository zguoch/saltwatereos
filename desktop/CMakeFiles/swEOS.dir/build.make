# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /Applications/CMake.app/Contents/bin/cmake

# The command to remove a file.
RM = /Applications/CMake.app/Contents/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop

# Include any dependencies generated for this target.
include CMakeFiles/swEOS.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/swEOS.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/swEOS.dir/flags.make

ui_MainWindow.h: MainWindow.ui
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating ui_MainWindow.h"
	/Users/zguo/Qt/5.9.9/clang_64/bin/uic -o /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/ui_MainWindow.h /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/MainWindow.ui

qrc_icons.cpp: Icons/fileopen.png
qrc_icons.cpp: Icons/filesave.png
qrc_icons.cpp: Icons/print.png
qrc_icons.cpp: Icons/help.png
qrc_icons.cpp: Icons/icons.qrc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating qrc_icons.cpp"
	/Users/zguo/Qt/5.9.9/clang_64/bin/rcc --name icons --output /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/qrc_icons.cpp /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/Icons/icons.qrc

CMakeFiles/swEOS.dir/main.cxx.o: CMakeFiles/swEOS.dir/flags.make
CMakeFiles/swEOS.dir/main.cxx.o: main.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/swEOS.dir/main.cxx.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/swEOS.dir/main.cxx.o -c /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/main.cxx

CMakeFiles/swEOS.dir/main.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/swEOS.dir/main.cxx.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/main.cxx > CMakeFiles/swEOS.dir/main.cxx.i

CMakeFiles/swEOS.dir/main.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/swEOS.dir/main.cxx.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/main.cxx -o CMakeFiles/swEOS.dir/main.cxx.s

CMakeFiles/swEOS.dir/main.cxx.o.requires:

.PHONY : CMakeFiles/swEOS.dir/main.cxx.o.requires

CMakeFiles/swEOS.dir/main.cxx.o.provides: CMakeFiles/swEOS.dir/main.cxx.o.requires
	$(MAKE) -f CMakeFiles/swEOS.dir/build.make CMakeFiles/swEOS.dir/main.cxx.o.provides.build
.PHONY : CMakeFiles/swEOS.dir/main.cxx.o.provides

CMakeFiles/swEOS.dir/main.cxx.o.provides.build: CMakeFiles/swEOS.dir/main.cxx.o


CMakeFiles/swEOS.dir/MainWindow.cxx.o: CMakeFiles/swEOS.dir/flags.make
CMakeFiles/swEOS.dir/MainWindow.cxx.o: MainWindow.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/swEOS.dir/MainWindow.cxx.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/swEOS.dir/MainWindow.cxx.o -c /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/MainWindow.cxx

CMakeFiles/swEOS.dir/MainWindow.cxx.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/swEOS.dir/MainWindow.cxx.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/MainWindow.cxx > CMakeFiles/swEOS.dir/MainWindow.cxx.i

CMakeFiles/swEOS.dir/MainWindow.cxx.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/swEOS.dir/MainWindow.cxx.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/MainWindow.cxx -o CMakeFiles/swEOS.dir/MainWindow.cxx.s

CMakeFiles/swEOS.dir/MainWindow.cxx.o.requires:

.PHONY : CMakeFiles/swEOS.dir/MainWindow.cxx.o.requires

CMakeFiles/swEOS.dir/MainWindow.cxx.o.provides: CMakeFiles/swEOS.dir/MainWindow.cxx.o.requires
	$(MAKE) -f CMakeFiles/swEOS.dir/build.make CMakeFiles/swEOS.dir/MainWindow.cxx.o.provides.build
.PHONY : CMakeFiles/swEOS.dir/MainWindow.cxx.o.provides

CMakeFiles/swEOS.dir/MainWindow.cxx.o.provides.build: CMakeFiles/swEOS.dir/MainWindow.cxx.o


CMakeFiles/swEOS.dir/qrc_icons.cpp.o: CMakeFiles/swEOS.dir/flags.make
CMakeFiles/swEOS.dir/qrc_icons.cpp.o: qrc_icons.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/swEOS.dir/qrc_icons.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/swEOS.dir/qrc_icons.cpp.o -c /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/qrc_icons.cpp

CMakeFiles/swEOS.dir/qrc_icons.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/swEOS.dir/qrc_icons.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/qrc_icons.cpp > CMakeFiles/swEOS.dir/qrc_icons.cpp.i

CMakeFiles/swEOS.dir/qrc_icons.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/swEOS.dir/qrc_icons.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/qrc_icons.cpp -o CMakeFiles/swEOS.dir/qrc_icons.cpp.s

CMakeFiles/swEOS.dir/qrc_icons.cpp.o.requires:

.PHONY : CMakeFiles/swEOS.dir/qrc_icons.cpp.o.requires

CMakeFiles/swEOS.dir/qrc_icons.cpp.o.provides: CMakeFiles/swEOS.dir/qrc_icons.cpp.o.requires
	$(MAKE) -f CMakeFiles/swEOS.dir/build.make CMakeFiles/swEOS.dir/qrc_icons.cpp.o.provides.build
.PHONY : CMakeFiles/swEOS.dir/qrc_icons.cpp.o.provides

CMakeFiles/swEOS.dir/qrc_icons.cpp.o.provides.build: CMakeFiles/swEOS.dir/qrc_icons.cpp.o


CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o: CMakeFiles/swEOS.dir/flags.make
CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o: swEOS_autogen/mocs_compilation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o -c /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/swEOS_autogen/mocs_compilation.cpp

CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/swEOS_autogen/mocs_compilation.cpp > CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.i

CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/swEOS_autogen/mocs_compilation.cpp -o CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.s

CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o.requires:

.PHONY : CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o.requires

CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o.provides: CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o.requires
	$(MAKE) -f CMakeFiles/swEOS.dir/build.make CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o.provides.build
.PHONY : CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o.provides

CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o.provides.build: CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o


# Object files for target swEOS
swEOS_OBJECTS = \
"CMakeFiles/swEOS.dir/main.cxx.o" \
"CMakeFiles/swEOS.dir/MainWindow.cxx.o" \
"CMakeFiles/swEOS.dir/qrc_icons.cpp.o" \
"CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o"

# External object files for target swEOS
swEOS_EXTERNAL_OBJECTS =

swEOS: CMakeFiles/swEOS.dir/main.cxx.o
swEOS: CMakeFiles/swEOS.dir/MainWindow.cxx.o
swEOS: CMakeFiles/swEOS.dir/qrc_icons.cpp.o
swEOS: CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o
swEOS: CMakeFiles/swEOS.dir/build.make
swEOS: /usr/local/lib/libvtkViewsQt-8.2.1.dylib
swEOS: /usr/local/lib/libvtkGUISupportQt-8.2.1.dylib
swEOS: /usr/local/lib/libvtkViewsInfovis-8.2.1.dylib
swEOS: /usr/local/lib/libvtkChartsCore-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersImaging-8.2.1.dylib
swEOS: /usr/local/lib/libvtkInfovisLayout-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingLabel-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingContextOpenGL2-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingGL2PSOpenGL2-8.2.1.dylib
swEOS: /usr/local/lib/libvtkgl2ps-8.2.1.dylib
swEOS: /usr/local/lib/libvtkViewsContext2D-8.2.1.dylib
swEOS: /Users/zguo/Qt/5.9.9/clang_64/lib/QtWidgets.framework/QtWidgets
swEOS: /Users/zguo/Qt/5.9.9/clang_64/lib/QtGui.framework/QtGui
swEOS: /Users/zguo/Qt/5.9.9/clang_64/lib/QtCore.framework/QtCore
swEOS: /usr/local/lib/libvtkInfovisCore-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingOpenGL2-8.2.1.dylib
swEOS: /usr/local/lib/libvtkglew-8.2.1.dylib
swEOS: /usr/local/lib/libvtkViewsCore-8.2.1.dylib
swEOS: /usr/local/lib/libvtkInteractionWidgets-8.2.1.dylib
swEOS: /usr/local/lib/libvtkInteractionStyle-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersExtraction-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersStatistics-8.2.1.dylib
swEOS: /usr/local/lib/libvtkImagingFourier-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersHybrid-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersModeling-8.2.1.dylib
swEOS: /usr/local/lib/libvtkImagingGeneral-8.2.1.dylib
swEOS: /usr/local/lib/libvtkImagingSources-8.2.1.dylib
swEOS: /usr/local/lib/libvtkImagingHybrid-8.2.1.dylib
swEOS: /usr/local/lib/libvtkIOImage-8.2.1.dylib
swEOS: /usr/local/lib/libvtkDICOMParser-8.2.1.dylib
swEOS: /usr/local/lib/libvtkmetaio-8.2.1.dylib
swEOS: /usr/local/lib/libvtkjpeg-8.2.1.dylib
swEOS: /usr/local/lib/libvtkpng-8.2.1.dylib
swEOS: /usr/local/lib/libvtktiff-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingAnnotation-8.2.1.dylib
swEOS: /usr/local/lib/libvtkImagingColor-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingVolume-8.2.1.dylib
swEOS: /usr/local/lib/libvtkImagingCore-8.2.1.dylib
swEOS: /usr/local/lib/libvtkIOXML-8.2.1.dylib
swEOS: /usr/local/lib/libvtkIOXMLParser-8.2.1.dylib
swEOS: /usr/local/lib/libvtkIOCore-8.2.1.dylib
swEOS: /usr/local/lib/libvtkdoubleconversion-8.2.1.dylib
swEOS: /usr/local/lib/libvtklz4-8.2.1.dylib
swEOS: /usr/local/lib/libvtklzma-8.2.1.dylib
swEOS: /usr/local/lib/libvtkexpat-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingContext2D-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingFreeType-8.2.1.dylib
swEOS: /usr/local/lib/libvtkRenderingCore-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersSources-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersGeneral-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonComputationalGeometry-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonColor-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersGeometry-8.2.1.dylib
swEOS: /usr/local/lib/libvtkFiltersCore-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonExecutionModel-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonDataModel-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonTransforms-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonMisc-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonMath-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonSystem-8.2.1.dylib
swEOS: /usr/local/lib/libvtkCommonCore-8.2.1.dylib
swEOS: /usr/local/lib/libvtksys-8.2.1.dylib
swEOS: /usr/local/lib/libvtkfreetype-8.2.1.dylib
swEOS: /usr/local/lib/libvtkzlib-8.2.1.dylib
swEOS: CMakeFiles/swEOS.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable swEOS"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/swEOS.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/swEOS.dir/build: swEOS

.PHONY : CMakeFiles/swEOS.dir/build

CMakeFiles/swEOS.dir/requires: CMakeFiles/swEOS.dir/main.cxx.o.requires
CMakeFiles/swEOS.dir/requires: CMakeFiles/swEOS.dir/MainWindow.cxx.o.requires
CMakeFiles/swEOS.dir/requires: CMakeFiles/swEOS.dir/qrc_icons.cpp.o.requires
CMakeFiles/swEOS.dir/requires: CMakeFiles/swEOS.dir/swEOS_autogen/mocs_compilation.cpp.o.requires

.PHONY : CMakeFiles/swEOS.dir/requires

CMakeFiles/swEOS.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/swEOS.dir/cmake_clean.cmake
.PHONY : CMakeFiles/swEOS.dir/clean

CMakeFiles/swEOS.dir/depend: ui_MainWindow.h
CMakeFiles/swEOS.dir/depend: qrc_icons.cpp
	cd /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop /Users/zguo/MyData/Research/3-CodeProject/EOS_H2ONaCl/desktop/CMakeFiles/swEOS.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/swEOS.dir/depend

