# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

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
CMAKE_COMMAND = /home/ignacio/clion-2019.3.2/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /home/ignacio/clion-2019.3.2/bin/cmake/linux/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/sequential_recalentamiento.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/sequential_recalentamiento.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/sequential_recalentamiento.dir/flags.make

CMakeFiles/sequential_recalentamiento.dir/main.cpp.o: CMakeFiles/sequential_recalentamiento.dir/flags.make
CMakeFiles/sequential_recalentamiento.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/sequential_recalentamiento.dir/main.cpp.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/sequential_recalentamiento.dir/main.cpp.o -c /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/main.cpp

CMakeFiles/sequential_recalentamiento.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/sequential_recalentamiento.dir/main.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/main.cpp > CMakeFiles/sequential_recalentamiento.dir/main.cpp.i

CMakeFiles/sequential_recalentamiento.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/sequential_recalentamiento.dir/main.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/main.cpp -o CMakeFiles/sequential_recalentamiento.dir/main.cpp.s

# Object files for target sequential_recalentamiento
sequential_recalentamiento_OBJECTS = \
"CMakeFiles/sequential_recalentamiento.dir/main.cpp.o"

# External object files for target sequential_recalentamiento
sequential_recalentamiento_EXTERNAL_OBJECTS =

sequential_recalentamiento: CMakeFiles/sequential_recalentamiento.dir/main.cpp.o
sequential_recalentamiento: CMakeFiles/sequential_recalentamiento.dir/build.make
sequential_recalentamiento: /usr/local/lib/libmongocxx.so
sequential_recalentamiento: /usr/local/lib/libbsoncxx.so
sequential_recalentamiento: CMakeFiles/sequential_recalentamiento.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable sequential_recalentamiento"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/sequential_recalentamiento.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/sequential_recalentamiento.dir/build: sequential_recalentamiento

.PHONY : CMakeFiles/sequential_recalentamiento.dir/build

CMakeFiles/sequential_recalentamiento.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/sequential_recalentamiento.dir/cmake_clean.cmake
.PHONY : CMakeFiles/sequential_recalentamiento.dir/clean

CMakeFiles/sequential_recalentamiento.dir/depend:
	cd /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug /home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/CMakeFiles/sequential_recalentamiento.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/sequential_recalentamiento.dir/depend
