# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.23

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = D:\software\path\CMake\bin\cmake.exe

# The command to remove a file.
RM = D:\software\path\CMake\bin\cmake.exe -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = D:\project\vessel

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = D:\project\vessel\build

# Include any dependencies generated for this target.
include CMakeFiles/LION.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/LION.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/LION.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/LION.dir/flags.make

CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.obj: CMakeFiles/LION.dir/flags.make
CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.obj: CMakeFiles/LION.dir/includes_CXX.rsp
CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.obj: LION_autogen/mocs_compilation.cpp
CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.obj: CMakeFiles/LION.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\project\vessel\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.obj"
	D:\software\path\MinGW\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.obj -MF CMakeFiles\LION.dir\LION_autogen\mocs_compilation.cpp.obj.d -o CMakeFiles\LION.dir\LION_autogen\mocs_compilation.cpp.obj -c D:\project\vessel\build\LION_autogen\mocs_compilation.cpp

CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.i"
	D:\software\path\MinGW\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\project\vessel\build\LION_autogen\mocs_compilation.cpp > CMakeFiles\LION.dir\LION_autogen\mocs_compilation.cpp.i

CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.s"
	D:\software\path\MinGW\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\project\vessel\build\LION_autogen\mocs_compilation.cpp -o CMakeFiles\LION.dir\LION_autogen\mocs_compilation.cpp.s

CMakeFiles/LION.dir/main.cpp.obj: CMakeFiles/LION.dir/flags.make
CMakeFiles/LION.dir/main.cpp.obj: CMakeFiles/LION.dir/includes_CXX.rsp
CMakeFiles/LION.dir/main.cpp.obj: ../main.cpp
CMakeFiles/LION.dir/main.cpp.obj: CMakeFiles/LION.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=D:\project\vessel\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/LION.dir/main.cpp.obj"
	D:\software\path\MinGW\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/LION.dir/main.cpp.obj -MF CMakeFiles\LION.dir\main.cpp.obj.d -o CMakeFiles\LION.dir\main.cpp.obj -c D:\project\vessel\main.cpp

CMakeFiles/LION.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/LION.dir/main.cpp.i"
	D:\software\path\MinGW\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E D:\project\vessel\main.cpp > CMakeFiles\LION.dir\main.cpp.i

CMakeFiles/LION.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/LION.dir/main.cpp.s"
	D:\software\path\MinGW\mingw64\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S D:\project\vessel\main.cpp -o CMakeFiles\LION.dir\main.cpp.s

# Object files for target LION
LION_OBJECTS = \
"CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.obj" \
"CMakeFiles/LION.dir/main.cpp.obj"

# External object files for target LION
LION_EXTERNAL_OBJECTS =

LION.exe: CMakeFiles/LION.dir/LION_autogen/mocs_compilation.cpp.obj
LION.exe: CMakeFiles/LION.dir/main.cpp.obj
LION.exe: CMakeFiles/LION.dir/build.make
LION.exe: libLION_UTILS.dll.a
LION.exe: libLION_GUI.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKLabelMap-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKFastMarching-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKPolynomials-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKBiasCorrection-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKColormap-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKConvolution-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKDICOMParser-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKDeformableMesh-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKDenoising-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKDiffusionTensorImage-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKPDEDeformableRegistration-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOBioRad-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOBruker-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOCSV-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOGE-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOHDF5-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOJPEG2000-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOLSM-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMINC-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMRC-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOSiemens-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOSpatialObjects-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOStimulate-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOTransformHDF5-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOTransformInsightLegacy-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOTransformMatlab-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKKLMRegionGrowing-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKMarkovRandomFieldsClassifiers-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKQuadEdgeMeshFiltering-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKRegionGrowing-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKRegistrationMethodsv4-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKTestKernel-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKVideoIO-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKVtkGlue-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKWatersheds-5.1.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkDomainsChemistryOpenGL2-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersFlowPaths-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersGeneric-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersHyperTree-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersParallelImaging-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersPoints-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersProgrammable-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersSMP-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersSelection-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersTopology-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersVerdict-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkGUISupportQtOpenGL-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkGUISupportQtSQL-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkGeovisCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOAMR-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOAsynchronous-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOCityGML-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOEnSight-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOExodus-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOExportOpenGL2-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOExportPDF-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOImport-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOInfovis-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOLSDyna-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOMINC-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOMovie-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOPLY-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOParallel-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOParallelXML-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOSegY-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOTecplotTable-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOVeraOut-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOVideo-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingMorphological-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingStatistics-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingStencil-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkInteractionImage-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingContextOpenGL2-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingImage-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingLOD-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingQt-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingVolumeOpenGL2-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkViewsContext2D-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkViewsQt-8.2.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_highgui455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_ml455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_objdetect455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_photo455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_stitching455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_video455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_videoio455.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKFFT-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkopenjpeg-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkminc2-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOIPL-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOXML-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/itkhdf5_cpp.lib
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/itkhdf5.lib
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOTransformBase-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKTransformFactory-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKImageFeature-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKOptimizersv4-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKOptimizers-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitklbfgs-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOBMP-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOGDCM-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkgdcmMSFF-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkgdcmDICT-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkgdcmIOD-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkgdcmDSED-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkgdcmCommon-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOGIPL-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOJPEG-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOTIFF-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitktiff-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkjpeg-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMeshBYU-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMeshFreeSurfer-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMeshGifti-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKgiftiio-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKEXPAT-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMeshOBJ-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMeshOFF-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMeshVTK-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMeshBase-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKQuadEdgeMesh-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOMeta-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKMetaIO-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIONIFTI-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKniftiio-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKznz-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIONRRD-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKNrrdIO-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOPNG-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkpng-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkzlib-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOVTK-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKIOImageBase-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKVideoCore-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKVTK-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKStatistics-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkNetlibSlatec-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKSpatialObjects-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKMesh-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKTransform-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKPath-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKCommon-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkdouble-conversion-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitksys-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKVNLInstantiation-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkvnl_algo-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkvnl-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkv3p_netlib-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libitkvcl-5.1.dll.a
LION.exe: D:/software/path/C++Package/ITK/itk-5.1.2-install-MinGW-x64-Release/lib/libITKSmoothing-5.1.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkDomainsChemistry-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkverdict-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOSQL-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtksqlite-8.2.dll.a
LION.exe: D:/software/path/Qt/Qt5.14.2/5.14.2/mingw73_64/lib/libQt5Sql.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkproj-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersAMR-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkpugixml-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOExport-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingGL2PSOpenGL2-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkgl2ps-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtklibharu-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtklibxml2-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtktheora-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkogg-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersParallel-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkexodusII-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOGeometry-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIONetCDF-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkNetCDF-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkjsoncpp-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkParallelCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOLegacy-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkhdf5_hl-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkhdf5-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersTexture-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingMath-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkGUISupportQt-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingOpenGL2-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkglew-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkViewsInfovis-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkChartsCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingContext2D-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersImaging-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkInfovisLayout-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkInfovisCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkViewsCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkInteractionWidgets-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersHybrid-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingGeneral-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingSources-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersModeling-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkInteractionStyle-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersExtraction-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersStatistics-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingFourier-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingHybrid-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOImage-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkDICOMParser-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkmetaio-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkpng-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtktiff-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkjpeg-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingAnnotation-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingColor-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingVolume-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkImagingCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOXML-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOXMLParser-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkIOCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkdoubleconversion-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtklz4-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtklzma-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkexpat-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingLabel-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingFreeType-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkRenderingCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonColor-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersGeometry-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersSources-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersGeneral-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonComputationalGeometry-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkFiltersCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonExecutionModel-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonDataModel-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonMisc-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonSystem-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonTransforms-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonMath-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkCommonCore-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtksys-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkfreetype-8.2.dll.a
LION.exe: D:/software/path/C++Package/VTK/VTK-8.2.0-install-MinGW-x64-Release/lib/libvtkzlib-8.2.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_calib3d455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_dnn455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_features2d455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_flann455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_imgcodecs455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_imgproc455.dll.a
LION.exe: D:/software/path/C++Package/OpenCV/opencv-4.5.5-install-MinGW-x64-Release/x64/mingw/lib/libopencv_core455.dll.a
LION.exe: D:/software/path/Qt/Qt5.14.2/5.14.2/mingw73_64/lib/libQt5Widgets.a
LION.exe: D:/software/path/Qt/Qt5.14.2/5.14.2/mingw73_64/lib/libQt5Gui.a
LION.exe: D:/software/path/Qt/Qt5.14.2/5.14.2/mingw73_64/lib/libQt5Core.a
LION.exe: CMakeFiles/LION.dir/linklibs.rsp
LION.exe: CMakeFiles/LION.dir/objects1.rsp
LION.exe: CMakeFiles/LION.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=D:\project\vessel\build\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable LION.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\LION.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/LION.dir/build: LION.exe
.PHONY : CMakeFiles/LION.dir/build

CMakeFiles/LION.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\LION.dir\cmake_clean.cmake
.PHONY : CMakeFiles/LION.dir/clean

CMakeFiles/LION.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" D:\project\vessel D:\project\vessel D:\project\vessel\build D:\project\vessel\build D:\project\vessel\build\CMakeFiles\LION.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/LION.dir/depend

