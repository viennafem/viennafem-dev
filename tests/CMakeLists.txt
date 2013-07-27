
#-----------------------------------------------------------------------------
# add the source files which should be tested without the trailing *.cpp
SET(PROGS poisson_2d poisson_3d integration)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
FOREACH(PROG ${PROGS})
   ADD_EXECUTABLE(${PROG}-test src/${PROG}.cpp)
   ADD_test(${PROG} ${PROG}-test)
ENDFOREACH(PROG)
#-----------------------------------------------------------------------------

ADD_EXECUTABLE(external_linkage  src/external_1.cpp
                                 src/external_2.cpp)
ADD_test(external_linkage external_linkage)
