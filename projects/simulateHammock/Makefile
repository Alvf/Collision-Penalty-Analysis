include ../include_top.mk

# calls:
INCLUDES   = ${INCLUDES_COMMON} -I../../ext/tinyobjloader
CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES}
LDFLAGS    = ${LDFLAGS_COMMON}

MAIN = simulateHammock.cpp
SOURCES    = ../../ext/glvu/glvu.cpp \
						 ../../src/util/MATRIX_UTIL.cpp \
						 ../../src/Collision/COLLISION_UTIL.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/DRAW_GL.cpp \
						 ../../src/util/PCG.cpp \
						 ../../src/util/DIAGONAL.cpp \
						 ../../src/util/BLOCK_SPARSE_MATRIX3.cpp \
						 ../../src/util/BLOCK_DIAGONAL_MATRIX3.cpp \
						 ../../src/Timestepper/Shell/TIMESTEPPER.cpp \
						 ../../src/Geometry/LINE_INTERSECT.cpp \
						 ../../src/Geometry/CUBE.cpp \
						 ../../src/Geometry/AABB_TREE.cpp \
						 ../../src/Geometry/TRIANGLE_MESH.cpp \
						 ../../src/Geometry/TRIANGLE_MESH_FASTER.cpp \
						 ../../src/Hyperelastic/Shell/THETA.cpp \
						 ../../src/Hyperelastic/Shell/THETA_FASTER.cpp \
						 ../../src/Hyperelastic/Shell/DIHEDRAL.cpp \
						 ../../src/Hyperelastic/Shell/ARAP.cpp \
						 ../../src/Hyperelastic/Shell/STRETCHING.cpp

include ../include_bottom.mk