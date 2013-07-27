# Project setup
###############

cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)


PROJECT(ViennaFEM)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

#specify custom boost path if needed:
#SET(BOOSTPATH "/opt/boost")
#SET(BOOSTPATH "C:/Program\ Files\ (x86)/boost/boost_1_42")

# alternate boost path for macs
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   INCLUDE_DIRECTORIES("/opt/local/include")
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

# specify include and source directory
################

INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES($ENV{VIENNAGRID})
INCLUDE_DIRECTORIES($ENV{VIENNAMATH})
INCLUDE_DIRECTORIES($ENV{VIENNADATA})
INCLUDE_DIRECTORIES($ENV{VIENNACL})
INCLUDE_DIRECTORIES(${BOOSTPATH})

# Common options
################

include(ViennaFEMCommon)

# Set release build:
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
ELSE()
 SET (CMAKE_BUILD_TYPE Release)
ENDIF()

# User options
##############

option(BUILD_DOXYGEN_DOCS "Build the Doxygen-generated API docs" OFF)

option(BUILD_MANUAL "Build the PDF manual" OFF)

cmake_dependent_option(ENABLE_DIST "Enable the the distribution tarball target" ON "BUILD_DOXYGEN_DOCS;BUILD_MANUAL" OFF)

mark_as_advanced(BUILD_DOXYGEN_DOCS BUILD_MANUAL ENABLE_DIST)

# Find prerequisites
####################

if(BUILD_DOXYGEN_DOCS)
   find_package(Doxygen REQUIRED)
endif()

if(BUILD_MANUAL)
   find_package(LATEX REQUIRED)
endif()

# Compile options
#################

include_directories(BEFORE ${PROJECT_SOURCE_DIR})


# Subdirectories
################

if(BUILD_EXAMPLES)
   add_subdirectory(examples)
endif()

if(BUILD_TESTING)
   INCLUDE(CTest)
   add_subdirectory(tests)
endif()

add_subdirectory(doc)

# Install
#########

install(DIRECTORY viennadata ${CMAKE_CURRENT_BINARY_DIR}/viennadata
   DESTINATION ${INSTALL_INCLUDE_DIR} COMPONENT dev
   FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")


# Add visibility of headers
# Necessary for Qt-Creator usage.
################
file( GLOB_RECURSE MAKE_HEADERS_VISIBLE_SRC
    *.hpp *.h)
add_custom_target( MAKE_HEADERS_VISIBLE SOURCES ${MAKE_HEADERS_VISIBLE_SRC} )



# packaging:
option(ENABLE_DIST "Enable the the distribution tarball target" OFF)
mark_as_advanced(ENABLE_DIST)

if(ENABLE_DIST)
   file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/dist")

   set(DISTDIR "${PROJECT_BINARY_DIR}/dist/ViennaFEM-${VERSION}")
   configure_file(CMakeLists.txt.in "${DISTDIR}/CMakeLists.txt" @ONLY)

   #
   # Full package:
   #
   add_custom_target(dist
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/LICENSE" "${DISTDIR}/LICENSE"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/README" "${DISTDIR}/README"

      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/build"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/build/README" "${DISTDIR}/build/README"

      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/doc"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/doc/viennacl.pdf" "${DISTDIR}/doc/viennacl.pdf"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/doc/viennadata.pdf" "${DISTDIR}/doc/viennadata.pdf"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/doc/viennafem.pdf" "${DISTDIR}/doc/viennafem.pdf"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/doc/viennagrid.pdf" "${DISTDIR}/doc/viennagrid.pdf"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/doc/viennamath.pdf" "${DISTDIR}/doc/viennamath.pdf"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/doc/doxygen"
         "${DISTDIR}/doc/doxygen"

      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/doc/Doxyfile" "${DISTDIR}/doc/Doxyfile"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/examples" "${DISTDIR}/examples"


      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/tests"
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/tests/src"
      COMMAND "${CMAKE_COMMAND}" -E make_directory "${DISTDIR}/tests/data"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/tests/src" "${DISTDIR}/tests/src"
      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/tests/data" "${DISTDIR}/tests/data"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/tests/clean.sh" "${DISTDIR}/tests/clean.sh"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/tests/CMakeLists.txt" "${DISTDIR}/tests/CMakeLists.txt"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/tests/CTestConfig.cmake" "${DISTDIR}/tests/CTestConfig.cmake"
      COMMAND "${CMAKE_COMMAND}" -E copy_if_different
         "${PROJECT_SOURCE_DIR}/tests/run_regression.sh" "${DISTDIR}/tests/run_regression.sh"



      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/viennacl" "${DISTDIR}/viennacl"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/viennadata" "${DISTDIR}/viennadata"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/viennafem" "${DISTDIR}/viennafem"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/viennagrid" "${DISTDIR}/viennagrid"

      COMMAND "${CMAKE_COMMAND}" -E copy_directory
         "${PROJECT_SOURCE_DIR}/viennamath" "${DISTDIR}/viennamath"


      COMMAND "${CMAKE_COMMAND}" -E tar czf ViennaFEM-${VERSION}.tar.gz
           ViennaFEM-${VERSION}
      WORKING_DIRECTORY "${PROJECT_BINARY_DIR}/dist"
      COMMENT "Generating dist/ViennaFEM-${VERSION}.tar.gz"
      VERBATIM)

endif()
