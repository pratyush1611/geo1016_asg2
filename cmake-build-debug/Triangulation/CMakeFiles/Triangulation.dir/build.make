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
include Triangulation\CMakeFiles\Triangulation.dir\depend.make

# Include the progress variables for this target.
include Triangulation\CMakeFiles\Triangulation.dir\progress.make

# Include the compile flags for this target's objects.
include Triangulation\CMakeFiles\Triangulation.dir\flags.make

Triangulation\CMakeFiles\Triangulation.dir\main.cpp.obj: Triangulation\CMakeFiles\Triangulation.dir\flags.make
Triangulation\CMakeFiles\Triangulation.dir\main.cpp.obj: ..\Triangulation\main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object Triangulation/CMakeFiles/Triangulation.dir/main.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\Triangulation.dir\main.cpp.obj /FdCMakeFiles\Triangulation.dir\ /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\main.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Triangulation.dir/main.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\Triangulation.dir\main.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\main.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Triangulation.dir/main.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\Triangulation.dir\main.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\main.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\triangulation.cpp.obj: Triangulation\CMakeFiles\Triangulation.dir\flags.make
Triangulation\CMakeFiles\Triangulation.dir\triangulation.cpp.obj: ..\Triangulation\triangulation.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object Triangulation/CMakeFiles/Triangulation.dir/triangulation.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\Triangulation.dir\triangulation.cpp.obj /FdCMakeFiles\Triangulation.dir\ /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\triangulation.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\triangulation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Triangulation.dir/triangulation.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\Triangulation.dir\triangulation.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\triangulation.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\triangulation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Triangulation.dir/triangulation.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\Triangulation.dir\triangulation.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\triangulation.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\triangulation_method.cpp.obj: Triangulation\CMakeFiles\Triangulation.dir\flags.make
Triangulation\CMakeFiles\Triangulation.dir\triangulation_method.cpp.obj: ..\Triangulation\triangulation_method.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object Triangulation/CMakeFiles/Triangulation.dir/triangulation_method.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\Triangulation.dir\triangulation_method.cpp.obj /FdCMakeFiles\Triangulation.dir\ /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\triangulation_method.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\triangulation_method.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Triangulation.dir/triangulation_method.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\Triangulation.dir\triangulation_method.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\triangulation_method.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\triangulation_method.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Triangulation.dir/triangulation_method.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\Triangulation.dir\triangulation_method.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\triangulation_method.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\matrix_algo.cpp.obj: Triangulation\CMakeFiles\Triangulation.dir\flags.make
Triangulation\CMakeFiles\Triangulation.dir\matrix_algo.cpp.obj: ..\Triangulation\matrix_algo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object Triangulation/CMakeFiles/Triangulation.dir/matrix_algo.cpp.obj"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoCMakeFiles\Triangulation.dir\matrix_algo.cpp.obj /FdCMakeFiles\Triangulation.dir\ /FS -c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\matrix_algo.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\matrix_algo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/Triangulation.dir/matrix_algo.cpp.i"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" > CMakeFiles\Triangulation.dir\matrix_algo.cpp.i @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\matrix_algo.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

Triangulation\CMakeFiles\Triangulation.dir\matrix_algo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/Triangulation.dir/matrix_algo.cpp.s"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\cl.exe" @<<
 /nologo /TP $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) /FoNUL /FAs /FaCMakeFiles\Triangulation.dir\matrix_algo.cpp.s /c E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation\matrix_algo.cpp
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

# Object files for target Triangulation
Triangulation_OBJECTS = \
"CMakeFiles\Triangulation.dir\main.cpp.obj" \
"CMakeFiles\Triangulation.dir\triangulation.cpp.obj" \
"CMakeFiles\Triangulation.dir\triangulation_method.cpp.obj" \
"CMakeFiles\Triangulation.dir\matrix_algo.cpp.obj"

# External object files for target Triangulation
Triangulation_EXTERNAL_OBJECTS =

bin\Triangulation.exe: Triangulation\CMakeFiles\Triangulation.dir\main.cpp.obj
bin\Triangulation.exe: Triangulation\CMakeFiles\Triangulation.dir\triangulation.cpp.obj
bin\Triangulation.exe: Triangulation\CMakeFiles\Triangulation.dir\triangulation_method.cpp.obj
bin\Triangulation.exe: Triangulation\CMakeFiles\Triangulation.dir\matrix_algo.cpp.obj
bin\Triangulation.exe: Triangulation\CMakeFiles\Triangulation.dir\build.make
bin\Triangulation.exe: lib\easy3d_core.lib
bin\Triangulation.exe: lib\easy3d_viewer.lib
bin\Triangulation.exe: lib\easy3d_optimizer.lib
bin\Triangulation.exe: lib\3rd_cminpack.lib
bin\Triangulation.exe: lib\3rd_glew.lib
bin\Triangulation.exe: lib\3rd_glfw.lib
bin\Triangulation.exe: lib\easy3d_fileio.lib
bin\Triangulation.exe: lib\easy3d_core.lib
bin\Triangulation.exe: lib\easy3d_util.lib
bin\Triangulation.exe: lib\3rd_glog.lib
bin\Triangulation.exe: Triangulation\CMakeFiles\Triangulation.dir\objects1.rsp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Linking CXX executable ..\bin\Triangulation.exe"
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	"D:\SW\clion\JetBrains\CLion 2020.3.2\bin\cmake\win\bin\cmake.exe" -E vs_link_exe --intdir=CMakeFiles\Triangulation.dir --rc="D:\Windows Kits\10\bin\10.0.18362.0\x86\rc.exe" --mt="D:\Windows Kits\10\bin\10.0.18362.0\x86\mt.exe" --manifests  -- "D:\SW\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.28.29333\bin\Hostx86\x86\link.exe" /nologo @CMakeFiles\Triangulation.dir\objects1.rsp @<<
 /out:..\bin\Triangulation.exe /implib:..\lib\Triangulation.lib /pdb:E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\bin\Triangulation.pdb /version:0.0  /machine:X86 /debug /INCREMENTAL /subsystem:console  ..\lib\easy3d_core.lib ..\lib\easy3d_viewer.lib ..\lib\easy3d_optimizer.lib ..\lib\3rd_cminpack.lib ..\lib\3rd_glew.lib opengl32.lib glu32.lib ..\lib\3rd_glfw.lib ..\lib\easy3d_fileio.lib ..\lib\easy3d_core.lib ..\lib\easy3d_util.lib ..\lib\3rd_glog.lib kernel32.lib user32.lib gdi32.lib winspool.lib shell32.lib ole32.lib oleaut32.lib uuid.lib comdlg32.lib advapi32.lib 
<<
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug

# Rule to build all files generated by this target.
Triangulation\CMakeFiles\Triangulation.dir\build: bin\Triangulation.exe

.PHONY : Triangulation\CMakeFiles\Triangulation.dir\build

Triangulation\CMakeFiles\Triangulation.dir\clean:
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation
	$(CMAKE_COMMAND) -P CMakeFiles\Triangulation.dir\cmake_clean.cmake
	cd E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug
.PHONY : Triangulation\CMakeFiles\Triangulation.dir\clean

Triangulation\CMakeFiles\Triangulation.dir\depend:
	$(CMAKE_COMMAND) -E cmake_depends "NMake Makefiles" E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\Triangulation E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation E:\TU_Delft\yr_1\q3\GEO1016\Assignment\asg02\A2_Triangulation_Code\cmake-build-debug\Triangulation\CMakeFiles\Triangulation.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : Triangulation\CMakeFiles\Triangulation.dir\depend

