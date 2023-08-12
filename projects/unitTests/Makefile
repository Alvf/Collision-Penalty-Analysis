include ../include_top.mk

# calls:
INCLUDES   = ${INCLUDES_COMMON} -I../../ext/catch2 -I../../ext/solvePoly
CFLAGS     = ${CFLAGS_COMMON} ${INCLUDES} -std=c++17
LDFLAGS    = ${LDFLAGS_COMMON} -std=c++17
MAIN			 = unitTests.cpp

SOURCES    = ../../ext/catch2/catch_amalgamated.cpp \
             ../../ext/solvePoly/poly34.cpp \
						 ../../src/Collision/COLLISION_UTIL.cpp \
						 ../../src/util/MATRIX_UTIL.cpp \
						 ../../src/util/MIN_HEAP.cpp \
						 ../../src/util/TIMER.cpp \
						 ../../src/util/FILE_IO.cpp \
						 ../../src/util/PCG.cpp \
						 ../../src/util/DIAGONAL.cpp \
						 ../../src/util/STRAND_DIAGONAL.cpp \
						 ../../src/util/BLOCK_SPARSE_MATRIX3.cpp \
						 ../../src/util/BLOCK_DIAGONAL_MATRIX3.cpp \
						 ../../src/Timestepper/Volume/QUASISTATIC.cpp \
						 ../../src/Timestepper/Volume/NEWMARK.cpp \
						 ../../src/Timestepper/Volume/BACKWARD_EULER_VELOCITY.cpp \
						 ../../src/Timestepper/Volume/BACKWARD_EULER_POSITION.cpp \
						 ../../src/Timestepper/Volume/TIMESTEPPER.cpp \
						 ../../src/Timestepper/Strand/TIMESTEPPER.cpp \
						 ../../src/Geometry/LINE_INTERSECT.cpp \
						 ../../src/Geometry/FIELD_3D.cpp \
						 ../../src/Geometry/FIELD_2D.cpp \
						 ../../src/Geometry/CAPSULE.cpp \
						 ../../src/Geometry/CYLINDER.cpp \
						 ../../src/Geometry/CUBE.cpp \
						 ../../src/Geometry/SPHERE.cpp \
						 ../../src/Geometry/STRAND_MESH.cpp \
						 ../../src/Geometry/STRAND_MESH_FASTER.cpp \
						 ../../src/Geometry/TET_MESH.cpp \
						 ../../src/Geometry/TET_MESH_FASTER.cpp \
						 ../../src/Geometry/AABB_TREE.cpp \
						 ../../src/Damping/Volume/ARAP_DAMPING.cpp \
						 ../../src/Damping/Volume/GREEN_DAMPING.cpp \
						 ../../src/Damping/Volume/EBERLE_DAMPING.cpp \
						 ../../src/Damping/Volume/EDGE_DAMPING.cpp \
						 ../../src/Hyperelastic/Volume/ANISOTROPIC_ARAP.cpp \
						 ../../src/Hyperelastic/Volume/ANISOTROPIC_STVK.cpp \
						 ../../src/Hyperelastic/Volume/ANISOTROPIC_FUNG.cpp \
						 ../../src/Hyperelastic/Volume/ANISOTROPIC_DIRICHLET.cpp \
						 ../../src/Hyperelastic/Volume/POWER_DIRICHLET.cpp \
						 ../../src/Hyperelastic/Volume/DIRICHLET.cpp \
						 ../../src/Hyperelastic/Volume/NEGATIVE_DIRICHLET.cpp \
						 ../../src/Hyperelastic/Volume/STVK.cpp \
						 ../../src/Hyperelastic/Volume/SNH_WITH_BARRIER.cpp \
						 ../../src/Hyperelastic/Volume/SNH.cpp \
						 ../../src/Hyperelastic/Volume/NEO_HOOKEAN_BW.cpp \
						 ../../src/Hyperelastic/Volume/ARAP.cpp \
						 ../../src/Hyperelastic/Volume/LINEAR.cpp \
						 ../../src/Hyperelastic/Volume/HYPERELASTIC.cpp \
						 ../../src/Hyperelastic/Volume/TET_STRAND_TWIST.cpp \
						 ../../src/Hyperelastic/Strand/QUADRATIC_STRETCHING.cpp \
						 ../../src/Hyperelastic/Strand/ISOTROPIC_THETA.cpp \
						 ../../src/Hyperelastic/Strand/QUADRATIC_BENDING.cpp \
						 ../../src/Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.cpp \
						 ../../src/Hyperelastic/Strand/TAN_BENDING.cpp \
						 ../../src/Hyperelastic/Shell/BENDING_SPRING.cpp \
						 ../../src/Hyperelastic/Shell/DIHEDRAL.cpp \
						 ../../src/Hyperelastic/Shell/THETA_FASTER.cpp \
						 ../../src/Hyperelastic/Shell/THETA.cpp \
						 ../../src/Hyperelastic/Shell/ARAP.cpp \
						 ../../src/Hyperelastic/Shell/STVK.cpp \
						 ../../src/Hyperelastic/Shell/BW_STRETCH.cpp \
						 ../../src/Hyperelastic/Shell/BW_SHEAR.cpp \
						 ../../src/Hyperelastic/Shell/BARAFF_WITKIN.cpp \
						 ../../src/Hyperelastic/Shell/STRETCHING.cpp

include ../include_bottom.mk
