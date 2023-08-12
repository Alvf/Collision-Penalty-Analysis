/*
This file is part of HOBAK.

HOBAK is free software: you can redistribute it and/or modify it under the terms of 
the GNU General Public License as published by the Free Software Foundation, either 
version 3 of the License, or (at your option) any later version.

HOBAK is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with HOBAK. 
If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef DRAW_GL_H
#define DRAW_GL_H

#include "Geometry/CAPSULE.h"
#include "Geometry/CYLINDER.h"
#include "Geometry/CUBE.h"
#include "Geometry/SPHERE.h"
#include "Geometry/AABB_TREE.h"
#include "Geometry/BOWL.h"
#include "Timestepper/Volume/TIMESTEPPER.h"
#include "Collision/ENERGY_12D.h"
#include "Collision/COLLISION_UTIL.h"
#include <glvu.h>

// DEBUG: here while TRIANGLE_MESH and SHELL::TIMESTEPPER are prototyping classes
#include "Geometry/TET_STRAND_MESH.h"
#include "Geometry/STRAND_MESH.h"
#include "Geometry/WISP_MESH.h"
#include "Geometry/TET_WISP_MESH.h"
#include "Geometry/TRIANGLE_MESH.h"
#include "Timestepper/Shell/TIMESTEPPER.h"
#include "Timestepper/Strand/TIMESTEPPER.h"

namespace HOBAK {

// Print a string to the GL window
void printGlString(string output);

void drawPlaneConstraints(const TET_MESH* tetMesh, const VOLUME::TIMESTEPPER* stepper);
void drawTet(const TET_MESH& mesh, const int tetID);
VECTOR3 planeNormal(const vector<VECTOR3>& plane);

// DEBUG: here while TRIANGLE_MESH and SHELL::TIMESTEPPER are prototyping classes
void drawKinematicConstraints(const TRIANGLE_MESH* triangleMesh, const SHELL::TIMESTEPPER* stepper);
void drawPlaneConstraints(const TRIANGLE_MESH* triangleMesh, const SHELL::TIMESTEPPER* stepper);
void drawKinematicConstraints(const STRAND_MESH* strandMesh, const STRAND::TIMESTEPPER* stepper);
void drawPlaneConstraints(const STRAND_MESH* strandMesh, const STRAND::TIMESTEPPER* stepper);

// draw the collision tets that have been built for this mesh
void drawCollisionTets(const TET_MESH& mesh);

// just draw a single vertex of a mesh
void drawVertex(const TET_MESH& mesh, const int index);
void drawVertex(const TRIANGLE_MESH& mesh, const int index);

// just draw the vertices of a tet mesh
void drawVertices(const TET_MESH& mesh, const vector<int>& toDraw);

// draw the collision cell around a specific surface triangle
void drawSurfaceFaceCollisionCell(const TET_MESH& mesh, const ENERGY_12D& collisionEnergy, const int triangleIndex);

// draw a specific surface triangles of a tet mesh
void drawSurfaceFace(const TET_MESH& mesh, const int index);

// draw only the surface triangles of a tet mesh
void drawSurfaceTriangles(const TET_MESH& mesh, bool drawOutlines);
void drawSurfaceTriangles(const TET_MESH& mesh, bool drawOutlines, const VECTOR3& tetColor, const VECTOR3& outlineColor);
void drawTriangleMesh(const TRIANGLE_MESH& mesh, bool drawOutlines);

// draw strands, or strand information
void drawWispMesh(const WISP_MESH& mesh, const int highlight = -1, const bool drawVertices = false);
void drawWispMesh(const TET_WISP_MESH& mesh, const int highlight = -1, const bool drawVertices = false);
void drawTetStrandMesh(const TET_STRAND_MESH& mesh, const int highlight = -1, const bool drawVertices = false);
void drawInvertedTetStrandMesh(const TET_STRAND_MESH& mesh);
void drawTet(const TET_STRAND_MESH& mesh, const int tetID);
void drawTetForces(const TET_STRAND_MESH& mesh, const int tetID);
void drawAnisotropy(const TET_STRAND_MESH& mesh, const int tetID);
void drawStrandMesh(const STRAND_MESH& mesh, const int highlight = -1, const bool drawVertices = false);
void drawStrandMeshOld(const STRAND_MESH& mesh, const int highlight = -1, const bool drawVertices = false);
void drawStrand(const STRAND_MESH& mesh, const int highlight = -1);
void drawCollisions(const STRAND_MESH& mesh);
void drawCollisionsOld(const STRAND_MESH& mesh);

// draw lines between differences between the two meshes
void drawDiff(const TET_MESH& mesh0, const TET_MESH& mesh1);

// draw kinematic objects
void drawKinematicShape(const KINEMATIC_SHAPE& shape);
void drawCapsule(const CAPSULE& capsule);
void drawCylinder(const CYLINDER& cylinder);
void drawSphere(const SPHERE& sphere);
void drawCube(const CUBE& cube);
void drawBowl(const BOWL& bowl);

// draw an AABB
void drawAABB(const VECTOR3& minCorner, const VECTOR3& maxCorner);
void drawAABB(const AABB_NODE& node);

// draw a AABB tree at a specific depth
void drawAABBTree(const AABB_NODE* node, const int drawDepth, const int currentDepth);
void drawAABBTree(const AABB_TREE& tree, const int drawDepth);

// draw coordinate axes, xyz = rgb
void drawAxes();

// draw collision information
void drawVertexFacePairs(const TET_MESH& tetMesh, const int highlighted = -1);
void drawEdgeEdgePairs(const TET_MESH& tetMesh, const int pairID = -1);

} // HOBAK

#endif
