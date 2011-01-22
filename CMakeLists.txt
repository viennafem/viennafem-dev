PROJECT(ViennaFEM)

cmake_minimum_required(VERSION 2.6)

#if you don't want the full compiler output, remove the following line
SET(CMAKE_VERBOSE_MAKEFILE ON)

#specify custom boost path if needed:
#SET(BOOSTPATH "/opt/boost")

# alternate boost path for macs
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
   INCLUDE_DIRECTORIES("/opt/local/include") 
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# ATI stream sdk:
# add include and lib path to build environment
SET(OPENCLPATH $ENV{ATISTREAMSDKROOT})
INCLUDE_DIRECTORIES(${OPENCLPATH}/include) 
LINK_DIRECTORIES(${OPENCLPATH}/lib/x86) 

#specify include and source directory
INCLUDE_DIRECTORIES(".")
INCLUDE_DIRECTORIES("../ViennaGrid/")
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
ADD_EXECUTABLE(tut1 examples/tutorial/tut1.cpp)

# Mac OS X specific linker part
#IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#   set(CMAKE_EXE_LINKER_FLAGS "-framework OpenCL")
#ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

# Linux specific linker part
#IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
#   TARGET_LINK_LIBRARIES(tut1 OpenCL)
#ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# Windows specific linker part
#IF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")
#   TARGET_LINK_LIBRARIES(tut1 OpenCL)
#ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Windows")

