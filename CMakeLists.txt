PROJECT(ViennaFEM)

set(VERSION_MAJOR 1)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)
set(VERSION ${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH})

cmake_minimum_required(VERSION 2.6)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

#specify custom boost path if needed:
#SET(BOOSTPATH "/opt/boost")
#SET(BOOSTPATH "C:/Program\ Files\ (x86)/boost/boost_1_42")

# alternate boost path for macs
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   INCLUDE_DIRECTORIES("/opt/local/include") 
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

#specify include and source directory
INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("../ViennaGrid/")
INCLUDE_DIRECTORIES("../ViennaData/")
INCLUDE_DIRECTORIES("../ViennaMath/")
INCLUDE_DIRECTORIES("../ViennaCL/")
INCLUDE_DIRECTORIES(${BOOSTPATH})

#add definitions, compiler switches, etc.
IF(DEFINED CMAKE_BUILD_TYPE)
 SET (CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE})
ELSE()
 SET (CMAKE_BUILD_TYPE Release)
ENDIF()

IF(CMAKE_COMPILER_IS_GNUCXX)
   #ADD_DEFINITIONS(-Wall -O0 -g -lpthread)
   ADD_DEFINITIONS(-Wall -pedantic -O3)
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

#list all source files here
ADD_EXECUTABLE(basic_poisson examples/tutorial/basic_poisson.cpp)
ADD_EXECUTABLE(poisson_1d examples/tutorial/poisson_1d.cpp)
#ADD_EXECUTABLE(poisson_2d examples/tutorial/poisson_2d.cpp)
#ADD_EXECUTABLE(poisson_2d_rect examples/tutorial/poisson_2d_rect.cpp)
#ADD_EXECUTABLE(poisson_cellquan_2d examples/tutorial/poisson_cellquan_2d.cpp)
#ADD_EXECUTABLE(poisson_3d examples/tutorial/poisson_3d.cpp)
#ADD_EXECUTABLE(poisson_3d_hex examples/tutorial/poisson_3d_hex.cpp)
#ADD_EXECUTABLE(lame_3d_hex examples/tutorial/lame_3d_hex.cpp)
#ADD_EXECUTABLE(lame_3d examples/tutorial/lame_3d.cpp)
#ADD_EXECUTABLE(sshape_2d examples/tutorial/sshape_2d.cpp)
#ADD_EXECUTABLE(sshape_3d examples/tutorial/sshape_3d.cpp)
#ADD_EXECUTABLE(weak_forms examples/tutorial/weak_forms.cpp)

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
