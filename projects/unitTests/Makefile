include ../include_top.mk

# calls:
INCLUDES   = ${INCLUDES_COMMON} -I../../ext/catch2 -I../../ext/solvePoly
CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES} -std=c++17
LDFLAGS    = ${LDFLAGS_COMMON} -std=c++17
MAIN			 = unitTests.cpp

SOURCES    = ../../ext/catch2/catch_amalgamated.cpp \
						 ../../src/Collision/COLLISION_UTIL.cpp \
						 ../../src/util/MATRIX_UTIL.cpp 

include ../include_bottom.mk
