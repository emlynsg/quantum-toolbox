# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.12

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
CMAKE_COMMAND = /home/exg103/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/182.3684.76/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/exg103/.local/share/JetBrains/Toolbox/apps/CLion/ch-0/182.3684.76/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = "/home/exg103/Documents/Honours/CPP Code"

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = "/home/exg103/Documents/Honours/CPP Code/cmake-build-release"

# Include any dependencies generated for this target.
include CMakeFiles/PracticeCode.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/PracticeCode.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/PracticeCode.dir/flags.make

CMakeFiles/PracticeCode.dir/main.cpp.o: CMakeFiles/PracticeCode.dir/flags.make
CMakeFiles/PracticeCode.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/exg103/Documents/Honours/CPP Code/cmake-build-release/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/PracticeCode.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PracticeCode.dir/main.cpp.o -c "/home/exg103/Documents/Honours/CPP Code/main.cpp"

CMakeFiles/PracticeCode.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PracticeCode.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/exg103/Documents/Honours/CPP Code/main.cpp" > CMakeFiles/PracticeCode.dir/main.cpp.i

CMakeFiles/PracticeCode.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PracticeCode.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/exg103/Documents/Honours/CPP Code/main.cpp" -o CMakeFiles/PracticeCode.dir/main.cpp.s

CMakeFiles/PracticeCode.dir/Grid.cpp.o: CMakeFiles/PracticeCode.dir/flags.make
CMakeFiles/PracticeCode.dir/Grid.cpp.o: ../Grid.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/exg103/Documents/Honours/CPP Code/cmake-build-release/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/PracticeCode.dir/Grid.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PracticeCode.dir/Grid.cpp.o -c "/home/exg103/Documents/Honours/CPP Code/Grid.cpp"

CMakeFiles/PracticeCode.dir/Grid.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PracticeCode.dir/Grid.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/exg103/Documents/Honours/CPP Code/Grid.cpp" > CMakeFiles/PracticeCode.dir/Grid.cpp.i

CMakeFiles/PracticeCode.dir/Grid.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PracticeCode.dir/Grid.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/exg103/Documents/Honours/CPP Code/Grid.cpp" -o CMakeFiles/PracticeCode.dir/Grid.cpp.s

CMakeFiles/PracticeCode.dir/Wavefunction.cpp.o: CMakeFiles/PracticeCode.dir/flags.make
CMakeFiles/PracticeCode.dir/Wavefunction.cpp.o: ../Wavefunction.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir="/home/exg103/Documents/Honours/CPP Code/cmake-build-release/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/PracticeCode.dir/Wavefunction.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/PracticeCode.dir/Wavefunction.cpp.o -c "/home/exg103/Documents/Honours/CPP Code/Wavefunction.cpp"

CMakeFiles/PracticeCode.dir/Wavefunction.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/PracticeCode.dir/Wavefunction.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E "/home/exg103/Documents/Honours/CPP Code/Wavefunction.cpp" > CMakeFiles/PracticeCode.dir/Wavefunction.cpp.i

CMakeFiles/PracticeCode.dir/Wavefunction.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/PracticeCode.dir/Wavefunction.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S "/home/exg103/Documents/Honours/CPP Code/Wavefunction.cpp" -o CMakeFiles/PracticeCode.dir/Wavefunction.cpp.s

# Object files for target PracticeCode
PracticeCode_OBJECTS = \
"CMakeFiles/PracticeCode.dir/main.cpp.o" \
"CMakeFiles/PracticeCode.dir/Grid.cpp.o" \
"CMakeFiles/PracticeCode.dir/Wavefunction.cpp.o"

# External object files for target PracticeCode
PracticeCode_EXTERNAL_OBJECTS =

PracticeCode: CMakeFiles/PracticeCode.dir/main.cpp.o
PracticeCode: CMakeFiles/PracticeCode.dir/Grid.cpp.o
PracticeCode: CMakeFiles/PracticeCode.dir/Wavefunction.cpp.o
PracticeCode: CMakeFiles/PracticeCode.dir/build.make
PracticeCode: CMakeFiles/PracticeCode.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir="/home/exg103/Documents/Honours/CPP Code/cmake-build-release/CMakeFiles" --progress-num=$(CMAKE_PROGRESS_4) "Linking CXX executable PracticeCode"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/PracticeCode.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/PracticeCode.dir/build: PracticeCode

.PHONY : CMakeFiles/PracticeCode.dir/build

CMakeFiles/PracticeCode.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/PracticeCode.dir/cmake_clean.cmake
.PHONY : CMakeFiles/PracticeCode.dir/clean

CMakeFiles/PracticeCode.dir/depend:
	cd "/home/exg103/Documents/Honours/CPP Code/cmake-build-release" && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" "/home/exg103/Documents/Honours/CPP Code" "/home/exg103/Documents/Honours/CPP Code" "/home/exg103/Documents/Honours/CPP Code/cmake-build-release" "/home/exg103/Documents/Honours/CPP Code/cmake-build-release" "/home/exg103/Documents/Honours/CPP Code/cmake-build-release/CMakeFiles/PracticeCode.dir/DependInfo.cmake" --color=$(COLOR)
.PHONY : CMakeFiles/PracticeCode.dir/depend

