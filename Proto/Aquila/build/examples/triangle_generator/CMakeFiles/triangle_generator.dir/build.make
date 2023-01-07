# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.23

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/quake/Projects/Current/SoundWave/c++/Aquila

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/quake/Projects/Current/SoundWave/c++/Aquila/build

# Include any dependencies generated for this target.
include examples/triangle_generator/CMakeFiles/triangle_generator.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include examples/triangle_generator/CMakeFiles/triangle_generator.dir/compiler_depend.make

# Include the progress variables for this target.
include examples/triangle_generator/CMakeFiles/triangle_generator.dir/progress.make

# Include the compile flags for this target's objects.
include examples/triangle_generator/CMakeFiles/triangle_generator.dir/flags.make

examples/triangle_generator/CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o: examples/triangle_generator/CMakeFiles/triangle_generator.dir/flags.make
examples/triangle_generator/CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o: /home/quake/Projects/Current/SoundWave/c++/Aquila/examples/triangle_generator/triangle_generator.cpp
examples/triangle_generator/CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o: examples/triangle_generator/CMakeFiles/triangle_generator.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/quake/Projects/Current/SoundWave/c++/Aquila/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/triangle_generator/CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o"
	cd /home/quake/Projects/Current/SoundWave/c++/Aquila/build/examples/triangle_generator && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT examples/triangle_generator/CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o -MF CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o.d -o CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o -c /home/quake/Projects/Current/SoundWave/c++/Aquila/examples/triangle_generator/triangle_generator.cpp

examples/triangle_generator/CMakeFiles/triangle_generator.dir/triangle_generator.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/triangle_generator.dir/triangle_generator.cpp.i"
	cd /home/quake/Projects/Current/SoundWave/c++/Aquila/build/examples/triangle_generator && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/quake/Projects/Current/SoundWave/c++/Aquila/examples/triangle_generator/triangle_generator.cpp > CMakeFiles/triangle_generator.dir/triangle_generator.cpp.i

examples/triangle_generator/CMakeFiles/triangle_generator.dir/triangle_generator.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/triangle_generator.dir/triangle_generator.cpp.s"
	cd /home/quake/Projects/Current/SoundWave/c++/Aquila/build/examples/triangle_generator && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/quake/Projects/Current/SoundWave/c++/Aquila/examples/triangle_generator/triangle_generator.cpp -o CMakeFiles/triangle_generator.dir/triangle_generator.cpp.s

# Object files for target triangle_generator
triangle_generator_OBJECTS = \
"CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o"

# External object files for target triangle_generator
triangle_generator_EXTERNAL_OBJECTS =

examples/triangle_generator/triangle_generator: examples/triangle_generator/CMakeFiles/triangle_generator.dir/triangle_generator.cpp.o
examples/triangle_generator/triangle_generator: examples/triangle_generator/CMakeFiles/triangle_generator.dir/build.make
examples/triangle_generator/triangle_generator: libAquila.a
examples/triangle_generator/triangle_generator: lib/libOoura_fft.a
examples/triangle_generator/triangle_generator: examples/triangle_generator/CMakeFiles/triangle_generator.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/quake/Projects/Current/SoundWave/c++/Aquila/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable triangle_generator"
	cd /home/quake/Projects/Current/SoundWave/c++/Aquila/build/examples/triangle_generator && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/triangle_generator.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/triangle_generator/CMakeFiles/triangle_generator.dir/build: examples/triangle_generator/triangle_generator
.PHONY : examples/triangle_generator/CMakeFiles/triangle_generator.dir/build

examples/triangle_generator/CMakeFiles/triangle_generator.dir/clean:
	cd /home/quake/Projects/Current/SoundWave/c++/Aquila/build/examples/triangle_generator && $(CMAKE_COMMAND) -P CMakeFiles/triangle_generator.dir/cmake_clean.cmake
.PHONY : examples/triangle_generator/CMakeFiles/triangle_generator.dir/clean

examples/triangle_generator/CMakeFiles/triangle_generator.dir/depend:
	cd /home/quake/Projects/Current/SoundWave/c++/Aquila/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/quake/Projects/Current/SoundWave/c++/Aquila /home/quake/Projects/Current/SoundWave/c++/Aquila/examples/triangle_generator /home/quake/Projects/Current/SoundWave/c++/Aquila/build /home/quake/Projects/Current/SoundWave/c++/Aquila/build/examples/triangle_generator /home/quake/Projects/Current/SoundWave/c++/Aquila/build/examples/triangle_generator/CMakeFiles/triangle_generator.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/triangle_generator/CMakeFiles/triangle_generator.dir/depend

