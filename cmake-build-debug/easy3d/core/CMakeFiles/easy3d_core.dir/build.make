# CMAKE generated file: DO NOT EDIT!
# Generated by "NMake Makefiles" Generator, CMake Version 3.17

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


.SUFFIXES: .hpux_make_needs_suffix_list


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE
NULL=nul
!ENDIF
SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "D:\SW\clion\JetBrains\CLion 2020.3.2\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "D:\SW\clion\JetBrains\CLion 2020.3.2\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

# Include any dependencies generated for this target.
include easy3d\core\CMakeFiles\easy3d_core.dir\depend.make

# Include the progress variables for this target.
include easy3d\core\CMakeFiles\easy3d_core.dir\progress.make

# Include the compile flags for this target's objects.
include easy3d\core\CMakeFiles\easy3d_core.dir\flags.make

easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.obj: ..\easy3d\core\graph.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/graph.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\graph.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\graph.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/graph.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\easy3d_core.dir\graph.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\graph.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/graph.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\graph.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\graph.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.obj: ..\easy3d\core\kdtree.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/kdtree.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\kdtree.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\kdtree.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/kdtree.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\easy3d_core.dir\kdtree.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\kdtree.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/kdtree.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\kdtree.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\kdtree.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.obj: ..\easy3d\core\point_cloud.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/point_cloud.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\point_cloud.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\point_cloud.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/point_cloud.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\easy3d_core.dir\point_cloud.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\point_cloud.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/point_cloud.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\point_cloud.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\point_cloud.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj: ..\easy3d\core\surface_mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/surface_mesh.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\surface_mesh.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/surface_mesh.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\easy3d_core.dir\surface_mesh.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\surface_mesh.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/surface_mesh.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\surface_mesh.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\surface_mesh.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj: easy3d\core\CMakeFiles\easy3d_core.dir\flags.make
easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj: ..\easy3d\core\manifold_builder.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object easy3d/core/CMakeFiles/easy3d_core.dir/manifold_builder.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj /FdCMakeFiles\easy3d_core.dir\easy3d_core.pdb /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\manifold_builder.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/easy3d_core.dir/manifold_builder.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\easy3d_core.dir\manifold_builder.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\manifold_builder.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/easy3d_core.dir/manifold_builder.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\easy3d_core.dir\manifold_builder.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core\manifold_builder.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

# Object files for target easy3d_core
easy3d_core_OBJECTS = \
"CMakeFiles\easy3d_core.dir\graph.cpp.obj" \
"CMakeFiles\easy3d_core.dir\kdtree.cpp.obj" \
"CMakeFiles\easy3d_core.dir\point_cloud.cpp.obj" \
"CMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj" \
"CMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj"

# External object files for target easy3d_core
easy3d_core_EXTERNAL_OBJECTS =

lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\graph.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\kdtree.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\point_cloud.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\surface_mesh.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\manifold_builder.cpp.obj
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\build.make
lib\easy3d_core.lib: easy3d\core\CMakeFiles\easy3d_core.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX static library ..\..\lib\easy3d_core.lib"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	$(CMAKE_COMMAND) -P CMakeFiles\easy3d_core.dir\cmake_clean_target.cmake
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\link.exe" /lib /nologo /machine:X86 /out:..\..\lib\easy3d_core.lib @CMakeFiles\easy3d_core.dir\objects1.rsp 
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

# Rule to build all files generated by this target.
easy3d\core\CMakeFiles\easy3d_core.dir\build: lib\easy3d_core.lib

.PHONY : easy3d\core\CMakeFiles\easy3d_core.dir\build

easy3d\core\CMakeFiles\easy3d_core.dir\clean:
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core
	$(CMAKE_COMMAND) -P CMakeFiles\easy3d_core.dir\cmake_clean.cmake
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug
.PHONY : easy3d\core\CMakeFiles\easy3d_core.dir\clean

easy3d\core\CMakeFiles\easy3d_core.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\easy3d\core E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\easy3d\core\CMakeFiles\easy3d_core.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : easy3d\core\CMakeFiles\easy3d_core.dir\depend

