cmake_minimum_required(VERSION 2.8.9)

project(LAScarSegmenter)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE "http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Modules/LAScarSegmenter")
set(EXTENSION_CATEGORY "Examples")
set(EXTENSION_STATUS "Beta")
set(EXTENSION_CONTRIBUTORS "LiangJia Zhu (UAB), Yi Gao (BWH), Josh Cates (Utah), Alan Morris (Utah), Danny Perry (Utah), Greg Gardner (Utah), Rob MacLeod (Utah), Allen Tannenbaum (UAB) ")
set(EXTENSION_DESCRIPTION "This is an module performing the LA scar segmentation")
set(EXTENSION_ICONURL "http://viewvc.slicer.org/viewvc.cgi/Slicer4/trunk/Extensions/Testing/LAScarSegmenter/LAScarSegmenter.png?revision=21745&view=co")
set(EXTENSION_SCREENSHOTURLS "http://wiki.slicer.org/slicerWiki/images/a/ab/Slicer-r19441-LAScarSegmenter-screenshot.png http://wiki.slicer.org/slicerWiki/images/1/1e/Slicer-r19441-LAScarSegmenter-screenshot-2.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(CLILAScarSegmenter)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})

