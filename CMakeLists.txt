cmake_minimum_required(VERSION 2.8.9)

project(LAScarSegmenter)

#-----------------------------------------------------------------------------
set(EXTENSION_HOMEPAGE " http://www.slicer.org/slicerWiki/index.php/Documentation/Nightly/Extensions/LAScarSegmenter")
set(EXTENSION_CATEGORY "Segmentation")
set(EXTENSION_STATUS "Beta")
set(EXTENSION_CONTRIBUTORS "Liangjia Zhu (SBU), Yi Gao (UAB), Josh Cates (Utah), Alan Morris (Utah), Danny Perry (Utah), Greg Gardner (Utah), Rob MacLeod (Utah), Allen Tannenbaum (SBU) ")
set(EXTENSION_DESCRIPTION "This is an module performing the LA scar segmentation")
set(EXTENSION_ICONURL "http://www.slicer.org/slicerWiki/images/5/56/LAScarSegmenter.png")
set(EXTENSION_SCREENSHOTURLS "http://www.slicer.org/slicerWiki/images/e/e5/LAScarSegmenterScreenShot.png")

#-----------------------------------------------------------------------------
find_package(Slicer REQUIRED)
include(${Slicer_USE_FILE})

#-----------------------------------------------------------------------------
add_subdirectory(CLILAScarSegmenter)

#-----------------------------------------------------------------------------
include(${Slicer_EXTENSION_CPACK})

