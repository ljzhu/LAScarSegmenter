#-----------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.8.9)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# See http://cmake.org/cmake/help/cmake-2-8-docs.html#section_Policies for details
#-----------------------------------------------------------------------------
if(POLICY CMP0017)
  cmake_policy(SET CMP0017 OLD)
endif()

#-----------------------------------------------------------------------------
# Module name
#   - Usually identique to the directory name)
#   - You shouldn't call 'project(...)'. See SlicerConfig.cmake for more details.
#-----------------------------------------------------------------------------
set(EXTENSION_NAME LAScarSegmenter)

#-----------------------------------------------------------------------------
# Prerequisites
#-----------------------------------------------------------------------------
if(NOT Slicer_SOURCE_DIR)
  find_package(Slicer REQUIRED)

  # Additional C/CXX flags - Should be defined before including Slicer_USE_FILE
  #set(ADDITIONAL_C_FLAGS "" CACHE STRING "Additional ${EXTENSION_NAME} C Flags")
  #set(ADDITIONAL_CXX_FLAGS "" CACHE STRING "Additional ${EXTENSION_NAME} CXX Flags")

  include(${Slicer_USE_FILE})

  set(EXTENSION_DEPENDS NA) # Specified as a space separated list or 'NA' if any
  set(EXTENSION_LICENSE_FILE ${Slicer_LICENSE_FILE})
  set(EXTENSION_README_FILE ${Slicer_README_FILE})
  set(EXTENSION_LICENSE_SHORT_DESCRIPTION ${Slicer_LICENSE_SHORT_DESCRIPTION})

  #include(SlicerEnableExtensionTesting)
endif()

# Configure XML description
configure_file(
  ${EXTENSION_NAME}.xml.in
  ${CMAKE_CURRENT_BINARY_DIR}/${EXTENSION_NAME}.xml
  @ONLY
  )

# Add module sub libraries
set(cli_module_export_directive "Q_SLICER_QTMODULES_LOADABLEEXTENSIONTEMPLATE_EXPORT")

# Additional includes
#set(cli_module_include_directories)

# Additional source files
#set(cli_module_SRCS)

# Additional Target libraries
set(cli_module_target_libraries ${ITK_LIBRARIES})


#-----------------------------------------------------------------------------
# Build
#-----------------------------------------------------------------------------
SEMMacroBuildCLI(
  NAME ${EXTENSION_NAME}
  #LOGO_HEADER ${Slicer_SOURCE_DIR}/Resources/NAMICLogo.h
  TARGET_LIBRARIES ${cli_module_target_libraries}
  ADDITIONAL_SRCS
    SFLS.h
    SFLSSegmentor3D.h
    SFLSSegmentor3D.txx
    utilities.hxx
    wallSegmenter3D.h
    wallSegmenter3D.txx
    localChanVeseWallSegmenter3D.h
    localChanVeseWallSegmenter3D.txx
  #INCLUDE_DIRECTORIES ${cli_module_include_directories}
  ADDITIONAL_SRCS ${cli_module_SRCS}
  #EXECUTABLE_ONLY
  )
