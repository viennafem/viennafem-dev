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

#
# ATI Stream SDK:
# add include and lib path to build environment
# The linker will most likely pick the correct link directory
# Otherwise, feel free to comment the wrong link directory for your machine.
#
SET(OPENCLPATH $ENV{ATISTREAMSDKROOT})
INCLUDE_DIRECTORIES(${OPENCLPATH}/include) 
LINK_DIRECTORIES(${OPENCLPATH}/lib/x86)       #use this on 32 bit systems
LINK_DIRECTORIES(${OPENCLPATH}/lib/x86_64)    #use this on 64 bit systems


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
ADD_EXECUTABLE(poisson_2d examples/tutorial/poisson_2d.cpp)
ADD_EXECUTABLE(poisson_3d examples/tutorial/poisson_3d.cpp)
ADD_EXECUTABLE(basic_lame examples/tutorial/basic_lame.cpp)
ADD_EXECUTABLE(lame_3d examples/tutorial/lame_3d.cpp)
ADD_EXECUTABLE(sshape_2d examples/tutorial/sshape_2d.cpp)
ADD_EXECUTABLE(sshape_3d examples/tutorial/sshape_3d.cpp)
#TARGET_LINK_LIBRARIES(poisson OpenCL)


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

