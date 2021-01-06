# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.18

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


# Produce verbose output by default.
VERBOSE = 1

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado/build

# Include any dependencies generated for this target.
include CMakeFiles/paralelizado.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/paralelizado.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/paralelizado.dir/flags.make

CMakeFiles/paralelizado.dir/main.cu.o: CMakeFiles/paralelizado.dir/flags.make
CMakeFiles/paralelizado.dir/main.cu.o: ../main.cu
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/paralelizado.dir/main.cu.o"
	/opt/cuda/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -x cu -c /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado/main.cu -o CMakeFiles/paralelizado.dir/main.cu.o

CMakeFiles/paralelizado.dir/main.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/paralelizado.dir/main.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/paralelizado.dir/main.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/paralelizado.dir/main.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target paralelizado
paralelizado_OBJECTS = \
"CMakeFiles/paralelizado.dir/main.cu.o"

# External object files for target paralelizado
paralelizado_EXTERNAL_OBJECTS =

paralelizado: CMakeFiles/paralelizado.dir/main.cu.o
paralelizado: CMakeFiles/paralelizado.dir/build.make
paralelizado: ../UTM.cpp
paralelizado: CMakeFiles/paralelizado.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CUDA executable paralelizado"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/paralelizado.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/paralelizado.dir/build: paralelizado

.PHONY : CMakeFiles/paralelizado.dir/build

CMakeFiles/paralelizado.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/paralelizado.dir/cmake_clean.cmake
.PHONY : CMakeFiles/paralelizado.dir/clean

CMakeFiles/paralelizado.dir/depend:
	cd /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado/build /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado/build /home/ignacio/Proyectos/S3E2_SAP/recocido_simulado_paralelizado/build/CMakeFiles/paralelizado.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/paralelizado.dir/depend

