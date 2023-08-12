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
#include "util/DRAW_GL.h"
#include <random>
#include "Hyperelastic/Volume/ANISOTROPIC_ARAP.h"

///////////////////////////////////////////////////////////////////////
// A bunch fo drawing routines for HOBAK objects, all in one place
//
// This is the only file with a GL dependency, so to remove it,
// just don't make any of the draw calls listed here.
///////////////////////////////////////////////////////////////////////

#include <glvu.h>
#if __APPLE__
#include <GL/glut.h>
#elif __linux__
#include <GL/glut.h>
#else
#include <GL/freeglut.h>
#include <GL/glu.h>
#endif

#include <iostream>

namespace HOBAK {

using namespace std;

///////////////////////////////////////////////////////////////////////
// Print a string to the GL window
///////////////////////////////////////////////////////////////////////
void printGlString(string output)
{
  glColor4f(10.0f, 0.0f, 0.0f, 10.0f);
  // Make ensuing transforms affect the projection matrix
  glMatrixMode(GL_PROJECTION);

  // set the projection matrix to an orthographic view
  glLoadIdentity();
  //glOrtho(-halfZoom, halfZoom, -halfZoom, halfZoom, -10, 10);
  glOrtho(0,0,0,0, -10, 10);

  // set the matric mode back to modelview
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  // must set color before setting raster position, otherwise it won't take
  glColor4f(10.0f, 0.0f, 0.0f, 10.0f);

  // normalized screen coordinates (-0.5, 0.5), due to the glLoadIdentity
  //glRasterPos3f(-halfZoom* 0.95, -halfZoom* 0.95, 0);
  //glRasterPos3f(-0.5, -0.5, 0);
  //glRasterPos3f(-0.5, -0.5, 0);
  glRasterPos3f(-0.95, -0.95, 0);

  glColor4f(10.0f, 0.0f, 0.0f, 10.0f);
  for (unsigned int x = 0; x < output.size(); x++)
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, output[x]);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawPlaneConstraints(const TET_MESH* tetMesh, const VOLUME::TIMESTEPPER* stepper)
{
  const vector<VECTOR3>& vertices = tetMesh->vertices();
  const vector<PLANE_CONSTRAINT>& constraints = stepper->planeConstraints();

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = constraints[x];
    const KINEMATIC_SHAPE* shape = constraint.shape;

    int index = constraint.vertexID;
    VECTOR3 vertex = vertices[index];
    const VECTOR3& localClosest = constraint.localClosestPoint;
    const VECTOR3& localNormal = constraint.localNormal;

    VECTOR3 closestPoint = shape->localVertexToWorld(localClosest); 

    glBegin(GL_POINTS);
      glColor4f(10.0, 0.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();

    glBegin(GL_POINTS);
      glColor4f(0.0, 0.0, 10.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // connect to closest point
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 10.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // draw the normal too
    VECTOR3 normal = shape->localNormalToWorld(localNormal);
    normal *= 0.1;
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 0.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
      glVertex3f(closestPoint[0] + normal[0], closestPoint[1] + normal[1], closestPoint[2] + normal[2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawPlaneConstraints(const TRIANGLE_MESH* triangleMesh, const SHELL::TIMESTEPPER* stepper)
{
  const vector<VECTOR3>& vertices = triangleMesh->vertices();
  const vector<PLANE_CONSTRAINT>& constraints = stepper->planeConstraints();

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = constraints[x];
    const KINEMATIC_SHAPE* shape = constraint.shape;

    int index = constraint.vertexID;
    VECTOR3 vertex = vertices[index];
    const VECTOR3& localClosest = constraint.localClosestPoint;
    const VECTOR3& localNormal = constraint.localNormal;

    VECTOR3 closestPoint = shape->localVertexToWorld(localClosest); 

    glBegin(GL_POINTS);
      glColor4f(10.0, 0.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();

    glBegin(GL_POINTS);
      glColor4f(0.0, 0.0, 10.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // connect to closest point
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 10.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // draw the normal too
    VECTOR3 normal = shape->localNormalToWorld(localNormal);
    normal *= 0.1;
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 0.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
      glVertex3f(closestPoint[0] + normal[0], closestPoint[1] + normal[1], closestPoint[2] + normal[2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawPlaneConstraints(const STRAND_MESH* strandMesh, const STRAND::TIMESTEPPER* stepper)
{
  //const vector<VECTOR3>& vertices = strandMesh->vertices();
  const vector<VECTOR3>& vertices = strandMesh->verticesOld();
  const vector<PLANE_CONSTRAINT>& constraints = stepper->planeConstraints();

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const PLANE_CONSTRAINT& constraint = constraints[x];
    const KINEMATIC_SHAPE* shape = constraint.shape;

    int index = constraint.vertexID;
    VECTOR3 vertex = vertices[index];
    const VECTOR3& localClosest = constraint.localClosestPoint;
    const VECTOR3& localNormal = constraint.localNormal;

    VECTOR3 closestPoint = shape->localVertexToWorld(localClosest); 

    glBegin(GL_POINTS);
      glColor4f(10.0, 0.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();

    glBegin(GL_POINTS);
      glColor4f(0.0, 0.0, 10.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // connect to closest point
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
    glEnd();

    // draw the normal too
    VECTOR3 normal = shape->localNormalToWorld(localNormal);
    normal *= 0.1;
    glBegin(GL_LINES);
      glColor4f(10.0, 10.0, 0.0, 1.0);
      glVertex3f(closestPoint[0], closestPoint[1], closestPoint[2]);
      glVertex3f(closestPoint[0] + normal[0], closestPoint[1] + normal[1], closestPoint[2] + normal[2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawKinematicConstraints(const TRIANGLE_MESH* triangleMesh, const SHELL::TIMESTEPPER* stepper)
{
  const vector<VECTOR3>& vertices = triangleMesh->vertices();
  const vector<KINEMATIC_CONSTRAINT>& constraints = stepper->kinematicConstraints();

  glPointSize(10.0);

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = constraints[x];
    int index = constraint.vertexID;

    const VECTOR3& vertex = vertices[index];
    glBegin(GL_POINTS);
      glColor4f(0.0, 10.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawKinematicConstraints(const STRAND_MESH* strandMesh, const STRAND::TIMESTEPPER* stepper)
{
  const vector<VECTOR3>& vertices = strandMesh->vertices();
  const vector<KINEMATIC_CONSTRAINT>& constraints = stepper->kinematicConstraints();

  glPointSize(10.0);

  for (unsigned int x = 0; x < constraints.size(); x++)
  {
    const KINEMATIC_CONSTRAINT& constraint = constraints[x];
    int index = constraint.vertexID;

    glDisable(GL_DEPTH_TEST);
    const VECTOR3& vertex = vertices[index];
    glBegin(GL_POINTS);
      glColor4f(0.0, 10.0, 0.0, 1.0);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
    glEnd();
    glEnable(GL_DEPTH_TEST);
  }
}

///////////////////////////////////////////////////////////////////////
// draw a single anisotropy direction on the mesh
///////////////////////////////////////////////////////////////////////
void drawAnisotropy(const TET_STRAND_MESH& mesh, const int tetID)
{
  const vector<VECTOR4I>& tets = mesh.tets();
  const vector<VECTOR3>& vertices = mesh.vertices();

  // draw the tet
  VECTOR3 v[] = {vertices[tets[tetID][0]],
                 vertices[tets[tetID][1]],
                 vertices[tets[tetID][2]],
                 vertices[tets[tetID][3]]};

  const VECTOR3 center = (v[0] + v[1] + v[2]) / 3.0;

  const MATRIX3& F = mesh.volumeFs()[tetID];
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  const MATRIX3 R = U * V.transpose();
  const MATRIX3 S = V * Sigma.asDiagonal() * V.transpose();

  const vector<STRAND::COMPOSITE>& materials = mesh.materials();
  const vector<VOLUME::HYPERELASTIC*>& volumeEnergies = materials[tetID].volumeEnergies();
  const VOLUME::ANISOTROPIC_ARAP* arap = (VOLUME::ANISOTROPIC_ARAP*)volumeEnergies[4];

  const REAL sign = (F.determinant() >= 0.0) ? 1.0 : -1.0;

  //const VECTOR3 aHat = S * arap->anisotropyDirection();
  //const VECTOR3 aHat = F * arap->anisotropyDirection();
  const VECTOR3 aHat = sign * (F * arap->anisotropyDirection());
  //const VECTOR3 aHat = R * arap->anisotropyDirection();
  const VECTOR3 tip = center + aHat;

  glDisable(GL_DEPTH_TEST);
  glColor4f(0.0, 0.0, 1.0, 1.0);

  glBegin(GL_POINTS);
    glVertex3f(center[0], center[1], center[2]);
    glVertex3f(tip[0], tip[1], tip[2]);
  glEnd();

  glBegin(GL_LINES);
    glVertex3f(center[0], center[1], center[2]);
    glVertex3f(tip[0], tip[1], tip[2]);
  glEnd();

  glEnable(GL_DEPTH_TEST);
}

///////////////////////////////////////////////////////////////////////
// draw a single tet in the mesh
///////////////////////////////////////////////////////////////////////
void drawTet(const TET_STRAND_MESH& mesh, const int tetID)
{
  const vector<VECTOR4I>& tets = mesh.tets();
  const vector<VECTOR3>& vertices = mesh.vertices();

  // draw the tet
  VECTOR3 v[] = {vertices[tets[tetID][0]],
                 vertices[tets[tetID][1]],
                 vertices[tets[tetID][2]],
                 vertices[tets[tetID][3]]};

  //glColor4f(1.0, 0.0, 0.0, 0.1);

  VECTOR3 n;
  glBegin(GL_TRIANGLES);
  /*
    n = (v[1] - v[0]).cross(v[3] - v[0]).normalized();
    glNormal3f(n[0], n[1], n[2]);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
    */
    
    n = (v[2] - v[0]).cross(v[1] - v[0]).normalized();
    glNormal3f(n[0], n[1], n[2]);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    
    /*
    n = (v[3] - v[0]).cross(v[2] - v[0]).normalized();
    glNormal3f(n[0], n[1], n[2]);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    
    n = (v[2] - v[1]).cross(v[3] - v[1]).normalized();
    glNormal3f(n[0], n[1], n[2]);
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
    */
  glEnd();

  glLineWidth(2.0);
  glColor4f(0, 0, 0, 1.0);
  glBegin(GL_LINES);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]);

    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]);

    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 

    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
  glEnd();
  /*
  glBegin(GL_POINTS);
    glColor4f(10.0, 0.0, 0.0, 10.0);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glColor4f(0.0, 0.0, 10.0, 10.0);
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
  glEnd();
  */
}

///////////////////////////////////////////////////////////////////////
// draw a single tet in the mesh
///////////////////////////////////////////////////////////////////////
void drawTet(const TET_MESH& mesh, const int tetID)
{
  const vector<VECTOR4I>& tets = mesh.tets();
  const vector<VECTOR3>& vertices = mesh.vertices();

  // draw the tet
  VECTOR3 v[] = {vertices[tets[tetID][0]],
                 vertices[tets[tetID][1]],
                 vertices[tets[tetID][2]],
                 vertices[tets[tetID][3]]};

  VECTOR3 n;
  glBegin(GL_TRIANGLES);
    n = (v[1] - v[0]).cross(v[3] - v[0]).normalized();
    glNormal3f(n[0], n[1], n[2]);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
    
    n = (v[2] - v[0]).cross(v[1] - v[0]).normalized();
    glNormal3f(n[0], n[1], n[2]);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    
    n = (v[3] - v[0]).cross(v[2] - v[0]).normalized();
    glNormal3f(n[0], n[1], n[2]);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    
    n = (v[2] - v[1]).cross(v[3] - v[1]).normalized();
    glNormal3f(n[0], n[1], n[2]);
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
  glEnd();

  glLineWidth(2.0);
  glColor4f(0, 0, 0, 1.0);
  glBegin(GL_LINES);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]);

    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]);

    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 

    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
  glEnd();
  /*
  glBegin(GL_POINTS);
    glColor4f(10.0, 0.0, 0.0, 10.0);
    glVertex3f(v[0][0], v[0][1], v[0][2]); 
    glColor4f(0.0, 0.0, 10.0, 10.0);
    glVertex3f(v[1][0], v[1][1], v[1][2]); 
    glVertex3f(v[2][0], v[2][1], v[2][2]); 
    glVertex3f(v[3][0], v[3][1], v[3][2]); 
  glEnd();
  */
}

///////////////////////////////////////////////////////////////////////
// draw a single tet in the mesh
///////////////////////////////////////////////////////////////////////
void drawTetForces(const TET_STRAND_MESH& mesh, const int tetID)
{
  const vector<VECTOR4I>& tets = mesh.tets();
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<VECTOR12>& forces = mesh.perTetForces();

  // draw the forces
  VECTOR3 v[] = {vertices[tets[tetID][0]],
                 vertices[tets[tetID][1]],
                 vertices[tets[tetID][2]],
                 vertices[tets[tetID][3]]};

  // endpoints
  vector<VECTOR3> tips(4);
  const VECTOR12& force = forces[tetID];
  for (int x = 0; x < 4; x++)
  {
    VECTOR3 y;
    y = force.segment<3>(3 * x);
    //tips[x] = v[x] + y.normalized();
    //tips[x] = v[x] + 0.01 * y;
    tips[x] = v[x] + 1e-4 * y;
  }

  glColor4f(0,0,1.0, 1.0);
  glBegin(GL_POINTS);
    for (int x = 0; x < 4; x++)
    {
      glVertex3f(v[x][0], v[x][1], v[x][2]);
      glVertex3f(tips[x][0], tips[x][1], tips[x][2]);
    }
  glEnd();

  glColor4f(0,1,0, 1.0);
  glBegin(GL_LINES);
    for (int x = 0; x < 4; x++)
    {
      glVertex3f(v[x][0], v[x][1], v[x][2]);
      glVertex3f(tips[x][0], tips[x][1], tips[x][2]);
    }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the collision tets that have been built for this mesh
///////////////////////////////////////////////////////////////////////
void drawCollisionTets(const TET_MESH& mesh)
{
  const vector<VECTOR4I>& vertexFaceCollisionTets = mesh.vertexFaceCollisionTets();
  const vector<VECTOR3>& vertices = mesh.vertices();

  // draw the tets
  for (unsigned int x = 0; x < vertexFaceCollisionTets.size(); x++)
  {
    VECTOR3 v[] = {vertices[vertexFaceCollisionTets[x][0]],
                   vertices[vertexFaceCollisionTets[x][1]],
                   vertices[vertexFaceCollisionTets[x][2]],
                   vertices[vertexFaceCollisionTets[x][3]]};

    VECTOR3 n;

    glColor4f(1.0, 0.0, 1.0, 1.0);
    glBegin(GL_TRIANGLES);
      n = (v[1] - v[0]).cross(v[3] - v[0]).normalized();
      glNormal3f(n[0], n[1], n[2]);
      glVertex3f(v[0][0], v[0][1], v[0][2]); 
      glVertex3f(v[1][0], v[1][1], v[1][2]); 
      glVertex3f(v[3][0], v[3][1], v[3][2]); 
      
      n = (v[2] - v[0]).cross(v[1] - v[0]).normalized();
      glNormal3f(n[0], n[1], n[2]);
      glVertex3f(v[0][0], v[0][1], v[0][2]); 
      glVertex3f(v[2][0], v[2][1], v[2][2]); 
      glVertex3f(v[1][0], v[1][1], v[1][2]); 
      
      n = (v[3] - v[0]).cross(v[2] - v[0]).normalized();
      glNormal3f(n[0], n[1], n[2]);
      glVertex3f(v[0][0], v[0][1], v[0][2]); 
      glVertex3f(v[3][0], v[3][1], v[3][2]); 
      glVertex3f(v[2][0], v[2][1], v[2][2]); 
      
      n = (v[2] - v[1]).cross(v[3] - v[1]).normalized();
      glNormal3f(n[0], n[1], n[2]);
      glVertex3f(v[1][0], v[1][1], v[1][2]); 
      glVertex3f(v[2][0], v[2][1], v[2][2]); 
      glVertex3f(v[3][0], v[3][1], v[3][2]); 
    glEnd();

    glBegin(GL_POINTS);
      glColor4f(10.0, 0.0, 0.0, 10.0);
      glVertex3f(v[0][0], v[0][1], v[0][2]); 
      glColor4f(0.0, 0.0, 10.0, 10.0);
      glVertex3f(v[1][0], v[1][1], v[1][2]); 
      glVertex3f(v[2][0], v[2][1], v[2][2]); 
      glVertex3f(v[3][0], v[3][1], v[3][2]); 
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
// just draw a single vertex of a tet mesh
///////////////////////////////////////////////////////////////////////
void drawVertex(const TET_MESH& mesh, const int index)
{
  const vector<VECTOR3>& vertices = mesh.vertices();

  if (index < 0) return;
  if (index >= (int)vertices.size()) return;

  glBegin(GL_POINTS);
    const VECTOR3& v = vertices[index];
    glVertex3dv(v.data());
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// just draw a single vertex of a triangle mesh
///////////////////////////////////////////////////////////////////////
void drawVertex(const TRIANGLE_MESH& mesh, const int index)
{
  const vector<VECTOR3>& vertices = mesh.vertices();

  if (index < 0) return;
  if (index >= (int)vertices.size()) return;

  glBegin(GL_POINTS);
    const VECTOR3& v = vertices[index];
    glVertex3dv(v.data());
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// just draw the vertices of a tet mesh
///////////////////////////////////////////////////////////////////////
void drawVertices(const TET_MESH& mesh, const vector<int>& toDraw)
{
  const vector<VECTOR3>& vertices = mesh.vertices();

  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < toDraw.size(); x++)
  {
    const VECTOR3& v = vertices[toDraw[x]];
    glVertex3dv(v.data());
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR3 planeNormal(const vector<VECTOR3>& plane)
{
  const VECTOR3 edge1 = plane[1] - plane[0];
  const VECTOR3 edge2 = plane[2] - plane[0];
  return edge1.cross(edge2).normalized();
}

///////////////////////////////////////////////////////////////////////
// draw the collision cell around a specific surface triangle
///////////////////////////////////////////////////////////////////////
void drawSurfaceFaceCollisionCell(const TET_MESH& mesh, const ENERGY_12D& collisionEnergy, const int triangleIndex)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<VECTOR3I>& surfaceTriangles = mesh.surfaceTriangles();
  const vector<VECTOR3I>& surfaceTriangleNeighbors = mesh.surfaceTriangleNeighbors();

  if (triangleIndex >= (int)surfaceTriangles.size())
    return;

  const VECTOR3I& tri = surfaceTriangles[triangleIndex];
  vector<VECTOR3> v;
  v.push_back(vertices[tri[0]]);
  v.push_back(vertices[tri[1]]);
  v.push_back(vertices[tri[2]]);
  VECTOR3 n = planeNormal(v);

  // get the normals of the three adjacent faces
  vector<VECTOR3> nNeighbors;
  const VECTOR3I& neighbors = surfaceTriangleNeighbors[triangleIndex];
  for (int x = 0; x < 3; x++)
  {
    assert(neighbors[x] != -1);
    const VECTOR3I& tNeighbor = surfaceTriangles[neighbors[x]];
    vector<VECTOR3> vNeighbor;
    vNeighbor.push_back(vertices[tNeighbor[0]]);
    vNeighbor.push_back(vertices[tNeighbor[1]]);
    vNeighbor.push_back(vertices[tNeighbor[2]]);
    VECTOR3 nNeighbor = planeNormal(vNeighbor);

    nNeighbors.push_back(nNeighbor);
  }

  // get the average normals along each edge
  vector<VECTOR3> edgeNormals(3);
  for (int x = 0; x < 3; x++)
  {
    const VECTOR3 ne = (nNeighbors[x] + n).normalized();
    edgeNormals[x] = ne;
  }

  // compute the vertices of the extruded face
  vector<VECTOR3> vExtruded(3);
  const REAL collisionEps = collisionEnergy.eps();
  for (int x = 0; x < 3; x++)
  {
    // get the direction
    const VECTOR3 e0 = edgeNormals[x];
    const VECTOR3 e1 = edgeNormals[(x + 1) % 3];
    const REAL t0 = e0.dot(n);
    const REAL t1 = e1.dot(n);

    const VECTOR3 v0 = v[x];
    const VECTOR3 v1 = v[(x + 1) % 3];
    const VECTOR3 v2 = v[(x + 2) % 3];
    const VECTOR3 v3 = v0 + e0 * collisionEps * (1.0 / t0);
    const VECTOR3 v5 = v1 + e1 * collisionEps * (1.0 / t1);

    const VECTOR3 n0 = ((v1 - v0).cross(v3 - v0)).normalized();
    const VECTOR3 n1 = ((v2 - v1).cross(v5 - v1)).normalized();

    const VECTOR direction = (n0.cross(n1)).normalized();
    const REAL t = direction.dot(n);
    vExtruded[(x + 1) % 3] = v1 + collisionEps * (1.0 / t) * direction;
  }

  // draw front and back
  glDisable(GL_CULL_FACE);

  // draw the extruded cap
  glColor4f(0.0, 0.0, 1.0, 0.25);
  glBegin(GL_TRIANGLES);
    glVertex3f(v[0][0], v[0][1], v[0][2]);
    glVertex3f(v[1][0], v[1][1], v[1][2]);
    glVertex3f(v[2][0], v[2][1], v[2][2]);
  glEnd();

  glLineWidth(1.0);

  // draw the side walls and outlines
  for (int x = 0; x < 3; x++)
  {
    VECTOR3 v0 = v[x];
    VECTOR3 v1 = v[(x + 1) % 3];
    
    VECTOR3 ve0 = vExtruded[x];
    VECTOR3 ve1 = vExtruded[(x + 1) % 3];

    // draw the faces along each extruded edge
    glColor4f(0.0, 0.0, 10.0, 0.5);
    glBegin(GL_TRIANGLES);
      glVertex3f(v0[0], v0[1], v0[2]);
      glVertex3f(v1[0], v1[1], v1[2]);
      glVertex3f(ve1[0], ve1[1], ve1[2]);

      glVertex3f(ve1[0], ve1[1], ve1[2]);
      glVertex3f(ve0[0], ve0[1], ve0[2]);
      glVertex3f(v0[0], v0[1], v0[2]);
    glEnd();
   
    // draw a yellow outline along the edge
    glColor4f(10.0, 10.0, 0.0, 10.0);
    glBegin(GL_LINES);
      glVertex3f(ve0[0], ve0[1], ve0[2]);
      glVertex3f(ve1[0], ve1[1], ve1[2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
// draw a specific surface triangles of a tet mesh
///////////////////////////////////////////////////////////////////////
void drawSurfaceFace(const TET_MESH& mesh, const int index)
{
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);

  const vector<VECTOR3I>& triangles = mesh.surfaceTriangles();

  // if it's not a valid face, give up gracefully
  if (index < 0) return;
  if (index > (int)triangles.size() - 1) return;

  const VECTOR3I& tri = triangles[index];
  glBegin(GL_TRIANGLES);
    glColor4f(0.0, 0.0, 10.0, 10.0);
    VECTOR3 v[3];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);

    // get the normal
    VECTOR3 edge1 = v[1] - v[0];
    VECTOR3 edge2 = v[2] - v[0];

    VECTOR3 normal = edge1.cross(edge2).normalized();
    glNormal3f(normal[0], normal[1], normal[2]);
      
    for (int y = 0; y < 3; y++)
      glVertex3f(v[y][0], v[y][1], v[y][2]);
  glEnd();

  glLineWidth(2.0);
  glColor4f(10, 0, 10, 1.0);
  glBegin(GL_LINE_STRIP);
    for (int y = 0; y < 4; y++)
    glVertex3f(v[y % 3][0], v[y % 3][1], v[y % 3][2]);
  glEnd();
}

// draw lines between differences between the two meshes
void drawDiff(const TET_MESH& mesh0, const TET_MESH& mesh1)
{
  const vector<VECTOR3>& vertices0 = mesh0.vertices();
  const vector<VECTOR3>& vertices1 = mesh1.vertices();

  assert(vertices0.size() == vertices1.size());

  glColor4f(10.0, 0.0, 10.0, 10.0);
  glBegin(GL_LINES);
  for (unsigned int x = 0; x < vertices0.size(); x++)
  {
    const VECTOR3 diff = vertices0[x] - vertices1[x];

    if (diff.norm() > 1e-6)
    {
      const VECTOR3 v0 = vertices0[x];
      const VECTOR3 v1 = vertices1[x];
      glVertex3f(v0[0], v0[1], v0[2]);
      glVertex3f(v1[0], v1[1], v1[2]);
    }
  }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw only the surface triangles of a tet mesh
///////////////////////////////////////////////////////////////////////
void drawSurfaceTriangles(const TET_MESH& mesh, bool drawOutlines)
{
  drawSurfaceTriangles(mesh, drawOutlines, 
                       VECTOR3(0.5, 0.5, 0.5), VECTOR3(0,0,0));
}

///////////////////////////////////////////////////////////////////////
// draw only the surface triangles of a tet mesh
///////////////////////////////////////////////////////////////////////
void drawSurfaceTriangles(const TET_MESH& mesh, bool drawOutlines, 
                          const VECTOR3& tetColor, const VECTOR3& outlineColor)
{
  const vector<VECTOR3I>& triangles = mesh.surfaceTriangles();
  VECTOR3 v[3];

  // draw front-facing triangles
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  //glColor4f(0.5, 0.5, 0.5, 1.0);
  glColor4f(tetColor[0], tetColor[1], tetColor[2], 1.0);
  glBegin(GL_TRIANGLES);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);

    // get the normal
    VECTOR3 edge1 = v[1] - v[0];
    VECTOR3 edge2 = v[2] - v[0];

    VECTOR3 normal = edge1.cross(edge2).normalized();
    glNormal3f(normal[0], normal[1], normal[2]);
      
    for (int y = 0; y < 3; y++)
      glVertex3f(v[y][0], v[y][1], v[y][2]);
  }
  glEnd();
 
  // draw back-facing a different color 
  glCullFace(GL_FRONT);
  glColor4f(1.0, 0.0, 1.0, 1.0);
  glBegin(GL_TRIANGLES);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);

    // get the normal
    VECTOR3 edge1 = v[1] - v[0];
    VECTOR3 edge2 = v[2] - v[0];

    VECTOR3 normal = edge1.cross(edge2).normalized();
    glNormal3f(normal[0], normal[1], normal[2]);
      
    for (int y = 0; y < 3; y++)
      glVertex3f(v[y][0], v[y][1], v[y][2]);
  }
  glEnd();

  // see if we're done
  if (!drawOutlines) return;

  glCullFace(GL_BACK);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);
    glLineWidth(2.0);
    //glColor4f(0, 0, 0, 1.0);
    glColor4f(outlineColor[0], outlineColor[1], outlineColor[1], 1.0);
    glBegin(GL_LINE_STRIP);
      for (int y = 0; y < 4; y++)
        glVertex3f(v[y % 3][0], v[y % 3][1], v[y % 3][2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
// draw the triangles of a shell
///////////////////////////////////////////////////////////////////////
void drawTriangleMesh(const TRIANGLE_MESH& mesh, bool drawOutlines)
{
  const vector<VECTOR3I>& triangles = mesh.triangles();
  VECTOR3 v[3];

  const VECTOR3 triangleColor(0.5, 0.5, 0.5);
  const VECTOR3 outlineColor(0, 0, 0);

  // draw front-facing triangles
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glColor4f(triangleColor[0], triangleColor[1], triangleColor[2], 1.0);
  glBegin(GL_TRIANGLES);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);

    // get the normal
    VECTOR3 edge1 = v[1] - v[0];
    VECTOR3 edge2 = v[2] - v[0];

    VECTOR3 normal = edge1.cross(edge2).normalized();
    glNormal3f(normal[0], normal[1], normal[2]);
      
    for (int y = 0; y < 3; y++)
      glVertex3f(v[y][0], v[y][1], v[y][2]);
  }
  glEnd();
 
  // draw back-facing a different color 
  glCullFace(GL_FRONT);
  glColor4f(1.0, 0.0, 1.0, 1.0);
  glBegin(GL_TRIANGLES);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);

    // get the normal
    VECTOR3 edge1 = v[1] - v[0];
    VECTOR3 edge2 = v[2] - v[0];

    VECTOR3 normal = edge1.cross(edge2).normalized();
    glNormal3f(normal[0], normal[1], normal[2]);
      
    for (int y = 0; y < 3; y++)
      glVertex3f(v[y][0], v[y][1], v[y][2]);
  }
  glEnd();

  // see if we're done
  if (!drawOutlines) return;

  glCullFace(GL_BACK);
  for (unsigned int x = 0; x < triangles.size(); x++)
  {
    const VECTOR3I& tri = triangles[x];
    for (int y = 0; y < 3; y++)
      v[y] = mesh.vertex(tri[y]);
    glLineWidth(2.0);
    //glColor4f(0, 0, 0, 1.0);
    glColor4f(outlineColor[0], outlineColor[1], outlineColor[1], 1.0);
    glBegin(GL_LINE_STRIP);
      for (int y = 0; y < 4; y++)
        glVertex3f(v[y % 3][0], v[y % 3][1], v[y % 3][2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
// draw a single strand
///////////////////////////////////////////////////////////////////////
void drawStrand(const STRAND_MESH& mesh, const int highlight)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();

  if (highlight >= strandIndices.size())
    return;

  // draw the strands
  glColor4f(0,0,0, 1.0);
  glLineWidth(5.0);
  glBegin(GL_LINES);
    for (unsigned int x = 0; x < strandIndices[highlight].size() - 1; x++)
    {
      const int i0 = strandIndices[highlight][x];
      const int i1 = strandIndices[highlight][x + 1];
      const VECTOR3& v0 = vertices[i0];
      const VECTOR3& v1 = vertices[i1];
      glVertex3f(v0[0], v0[1], v0[2]);
      glVertex3f(v1[0], v1[1], v1[2]);
    }
  glEnd();

  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    glColor4f(1.0, 0.0, 0.0, 1.0);
    glNormal3f(1.0, 0,0);
    for (unsigned int x = 0; x < strandIndices[highlight].size(); x++)
    {
      const int i0 = strandIndices[highlight][x];
      const VECTOR3& v0 = vertices[i0];
      glVertex3f(v0[0], v0[1], v0[2]);
    }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the lines in a strand mesh
///////////////////////////////////////////////////////////////////////
void drawStrandMesh(const STRAND_MESH& mesh, const int highlight, const bool drawVertices)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();

  // draw the strands
  glColor4f(0,0,0, 1.0);
  // glLineWidth(5.0);
  glLineWidth(2.0);
  glBegin(GL_LINES);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
      for (unsigned int x = 0; x < strandIndices[y].size() - 1; x++)
      {
        const int i0 = strandIndices[y][x];
        const int i1 = strandIndices[y][x + 1];
        const VECTOR3& v0 = vertices[i0];
        const VECTOR3& v1 = vertices[i1];
        glVertex3f(v0[0], v0[1], v0[2]);
        glVertex3f(v1[0], v1[1], v1[2]);
      }
  glEnd();

  if (!drawVertices) return;

  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      if (y == highlight)
        glColor4f(1.0, 1.0, 0.0, 1.0);
      else
        glColor4f(1.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the tets of a strand mesh
///////////////////////////////////////////////////////////////////////
void drawTetStrandMesh(const TET_STRAND_MESH& mesh, const int highlight, const bool drawVertices)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();

  for (int x = 0; x < mesh.totalTets(); x++)
    drawTet(mesh, x);

  if (!drawVertices) return;

  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      if (y == highlight)
        glColor4f(1.0, 1.0, 0.0, 1.0);
      else
        glColor4f(1.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the inverted tets of a strand mesh
///////////////////////////////////////////////////////////////////////
void drawInvertedTetStrandMesh(const TET_STRAND_MESH& mesh)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();
  const vector<bool>& isTetInverted = mesh.isTetInverted();

  glColor4f(100.0, 0, 0, 1.0);
  for (int x = 0; x < mesh.totalTets(); x++)
    if (isTetInverted[x])
      drawTet(mesh, x);

  /*
  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      glColor4f(1.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
  */
}

///////////////////////////////////////////////////////////////////////
// draw the lines in a strand mesh
///////////////////////////////////////////////////////////////////////
void drawWispMesh(const TET_WISP_MESH& mesh, const int highlight, const bool drawVertices)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();
  const VECTOR& vertexRadii = mesh.vertexRadii();

  glColor4f(0.25,0,0, 1.0);
  GLUquadricObj* quadric;
  quadric = gluNewQuadric();
  for (unsigned int y = 0; y < strandIndices.size(); y++)
    for (unsigned int x = 0; x < strandIndices[y].size() - 1; x++)
    {
      const int i0 = strandIndices[y][x];
      const int i1 = strandIndices[y][x + 1];
      const VECTOR3& v0 = vertices[i0];
      const VECTOR3& v1 = vertices[i1];
      const VECTOR3& t = v0;
      //const Eigen::AngleAxis<GLfloat> R{ cylinder.rotation().cast<GLfloat>() };

      const REAL radius0 = vertexRadii[i0];
      const REAL radius1 = vertexRadii[i1];
      const REAL height = (v0 - v1).norm();
      const VECTOR3 zAxis(0,0,1);
      const VECTOR3 want = (v1 - v0).normalized();

      // here's the axis to rotate about
      const VECTOR3 axis = zAxis.cross(want);

      // here's the angle
      const REAL angle = acos(zAxis.dot(want));

      glPushMatrix();
        // apply the transform
        glTranslatef(t[0], t[1], t[2]);
        glRotatef((180.0/M_PI) * angle, axis[0], axis[1], axis[2]);

        // draw it along y axis instead of z axis
        //glRotatef(90.0, 1.0, 0.0, 0.0);
        //glTranslatef(0,0,-0.5 * height);

        // draw the end caps
        glPushMatrix();
          glTranslatef(0.0, 0.0, height);
          gluDisk(quadric, 0.0, radius0, 10, 2);
        glPopMatrix();
        glPushMatrix();
          glRotatef(180.0, 0,1.0,0);
          gluDisk(quadric, 0.0, radius1, 10, 2);
        glPopMatrix();

        // draw the cylinder wall
        gluCylinder(quadric, radius0, radius1, height, 10, 10);
      glPopMatrix();
    }
  gluDeleteQuadric(quadric);

  if (!drawVertices) return;

  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      if (y == highlight)
        glColor4f(1.0, 1.0, 0.0, 1.0);
      else
        glColor4f(0.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
}
/*
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();
  const VECTOR& vertexRadii = mesh.vertexRadii();

  glColor4f(0.25,0,0, 1.0);
  GLUquadricObj* quadric;
  quadric = gluNewQuadric();
  for (unsigned int y = 0; y < strandIndices.size(); y++)
    for (unsigned int x = 0; x < strandIndices[y].size() - 1; x++)
    {
      const int i0 = strandIndices[y][x];
      const int i1 = strandIndices[y][x + 1];
      const VECTOR3& v0 = vertices[i0];
      const VECTOR3& v1 = vertices[i1];
      const VECTOR3& t = v0;
      //const Eigen::AngleAxis<GLfloat> R{ cylinder.rotation().cast<GLfloat>() };

      const REAL radius0 = vertexRadii[i0];
      const REAL radius1 = vertexRadii[i1];
      const REAL height = (v0 - v1).norm();
      const VECTOR3 zAxis(0,0,1);
      const VECTOR3 want = (v1 - v0).normalized();

      // here's the axis to rotate about
      const VECTOR3 axis = zAxis.cross(want);

      // here's the angle
      const REAL angle = acos(zAxis.dot(want));
      const Eigen::AngleAxis<GLfloat> R{ rotationMatrix(axis, angle).cast<GLfloat>() };

      glPushMatrix();
        // apply the transform
        glTranslatef(t[0], t[1], t[2]);
        glRotatef((180.0/M_PI) * R.angle(), R.axis().x(), R.axis().y(), R.axis().z());

        // draw it along y axis instead of z axis
        //glRotatef(90.0, 1.0, 0.0, 0.0);
        //glTranslatef(0,0,-0.5 * height);

        // draw the end caps
        glPushMatrix();
          glTranslatef(0.0, 0.0, height);
          gluDisk(quadric, 0.0, radius0, 10, 2);
        glPopMatrix();
        glPushMatrix();
          glRotatef(180.0, 0,1.0,0);
          gluDisk(quadric, 0.0, radius1, 10, 2);
        glPopMatrix();

        // draw the cylinder wall
        gluCylinder(quadric, radius0, radius1, height, 10, 10);
      glPopMatrix();
    }
  gluDeleteQuadric(quadric);

  if (!drawVertices) return;

  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      if (y == highlight)
        glColor4f(1.0, 1.0, 0.0, 1.0);
      else
        glColor4f(0.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
}
*/

///////////////////////////////////////////////////////////////////////
// draw the lines in a strand mesh
///////////////////////////////////////////////////////////////////////
void drawWispMesh(const WISP_MESH& mesh, const int highlight, const bool drawVertices)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();
  const VECTOR& wispRadii = mesh.wispRadii();

  glColor4f(0.25,0,0, 1.0);
  GLUquadricObj* quadric;
  quadric = gluNewQuadric();
  for (unsigned int y = 0; y < strandIndices.size(); y++)
    for (unsigned int x = 0; x < strandIndices[y].size() - 1; x++)
    {
      const int i0 = strandIndices[y][x];
      const int i1 = strandIndices[y][x + 1];
      const VECTOR3& v0 = vertices[i0];
      const VECTOR3& v1 = vertices[i1];
      const VECTOR3& t = v0;
      //const Eigen::AngleAxis<GLfloat> R{ cylinder.rotation().cast<GLfloat>() };

      const REAL radius0 = wispRadii[i0];
      const REAL radius1 = wispRadii[i1];
      const REAL height = (v0 - v1).norm();
      const VECTOR3 zAxis(0,0,1);
      const VECTOR3 want = (v1 - v0).normalized();

      // here's the axis to rotate about
      const VECTOR3 axis = zAxis.cross(want);

      // here's the angle
      const REAL angle = acos(zAxis.dot(want));
      const Eigen::AngleAxis<GLfloat> R{ rotationMatrix(axis, angle).cast<GLfloat>() };

      glPushMatrix();
        // apply the transform
        glTranslatef(t[0], t[1], t[2]);
        glRotatef((180.0/M_PI) * R.angle(), R.axis().x(), R.axis().y(), R.axis().z());

        // draw it along y axis instead of z axis
        //glRotatef(90.0, 1.0, 0.0, 0.0);
        //glTranslatef(0,0,-0.5 * height);

        // draw the end caps
        glPushMatrix();
          glTranslatef(0.0, 0.0, height);
          gluDisk(quadric, 0.0, radius0, 10, 2);
        glPopMatrix();
        glPushMatrix();
          glRotatef(180.0, 0,1.0,0);
          gluDisk(quadric, 0.0, radius1, 10, 2);
        glPopMatrix();

        // draw the cylinder wall
        gluCylinder(quadric, radius0, radius1, height, 10, 10);
      glPopMatrix();
    }
  gluDeleteQuadric(quadric);

  if (!drawVertices) return;

  // draw the vertices
  glPointSize(10.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      if (y == highlight)
        glColor4f(1.0, 1.0, 0.0, 1.0);
      else
        glColor4f(0.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw the lines in a strand mesh
///////////////////////////////////////////////////////////////////////
void drawStrandMeshOld(const STRAND_MESH& mesh, const int highlight, const bool drawVertices)
{
  const vector<VECTOR3>& vertices = mesh.verticesOld();
  const vector<vector<int> >& strandIndices = mesh.strandIndices();

  std::mt19937 gen(123456);
  //std::uniform_real_distribution<double> twister(0.0, 1.0);
  std::uniform_real_distribution<double> twister(0.0, 0.25);

  // draw the strands
  glLineWidth(1.0);
  glBegin(GL_LINES);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      // color each strand randomly
      REAL r = twister(gen);
      REAL g = twister(gen);
      REAL b = twister(gen);
      //glColor4f(r,g,b, 0.5);
      glColor4f(r,r,r,1.0);
      for (unsigned int x = 0; x < strandIndices[y].size() - 1; x++)
      {
        const int i0 = strandIndices[y][x];
        const int i1 = strandIndices[y][x + 1];
        const VECTOR3& v0 = vertices[i0];
        const VECTOR3& v1 = vertices[i1];
        glVertex3f(v0[0], v0[1], v0[2]);
        glVertex3f(v1[0], v1[1], v1[2]);
      }
    }
  glEnd();

  if (!drawVertices) return;
  
  //glPointSize(10.0);
  glPointSize(3.0);
  glBegin(GL_POINTS);
    for (unsigned int y = 0; y < strandIndices.size(); y++)
    {
      glNormal3f(1.0, 0,0);
      if (y == highlight)
        glColor4f(1.0, 1.0, 0.0, 1.0);
      else
        glColor4f(1.0, 0.0, 0.0, 1.0);
      for (unsigned int x = 0; x < strandIndices[y].size(); x++)
      {
        const int i0 = strandIndices[y][x];
        const VECTOR3& v0 = vertices[i0];
        glVertex3f(v0[0], v0[1], v0[2]);
      }
    }
  glEnd();
}

/*
///////////////////////////////////////////////////////////////////////
// draw the twisting forces on a strand mesh
///////////////////////////////////////////////////////////////////////
void drawStrandTwistForces(const STRAND_MESH& mesh)
{
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<VECTOR2I>& edges = mesh.edges();
  const VECTOR twistForces = mesh.computeTwistingForces();

  glColor4f(0,1.0,0,1.0);
  glDisable(GL_DEPTH_TEST);
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    const VECTOR3 vertex = vertices[x];

    VECTOR3 force;
    for (unsigned int y = 0; y < 3; y++)
      force[y] = twistForces[3 * x + y];

    const VECTOR3 forceLine = vertex + force;

    glBegin(GL_LINES);
      glVertex3f(vertex[0], vertex[1], vertex[2]);
      glVertex3f(forceLine[0], forceLine[1], forceLine[2]);
    glEnd();
  }
  glEnable(GL_DEPTH_TEST);
}
*/

///////////////////////////////////////////////////////////////////////
// draw a capsule
///////////////////////////////////////////////////////////////////////
void drawCapsule(const CAPSULE& capsule)
{
  const VECTOR3& t = capsule.translation();
  const Eigen::AngleAxis<GLfloat> R{ capsule.rotation().cast<GLfloat>() };

  const REAL radius = capsule.radius();
  const REAL height = capsule.height();

  GLUquadricObj* quadric;
  quadric = gluNewQuadric();

  glPushMatrix();
    // apply the transform
    glTranslatef(t[0], t[1], t[2]);
    glRotatef((180.0/M_PI) * R.angle(), R.axis().x(), R.axis().y(), R.axis().z());

    // draw it along y axis instead of z axis
    glRotatef(90.0, 1.0, 0.0, 0.0);
    glTranslatef(0,0,-0.5 * height);

    // draw the end caps
    glPushMatrix();
      glTranslatef(0.0, 0.0, height);
      //gluDisk(quadric, 0.0, radius, 20, 2);
      glutSolidSphere(radius, 20, 20);
    glPopMatrix();
    glPushMatrix();
      //glRotatef(180.0, 0,1.0,0);
      //gluDisk(quadric, 0.0, radius, 20, 2);
      glutSolidSphere(radius, 20, 20);
    glPopMatrix();

    // draw the cylinder wall
    gluCylinder(quadric, radius, radius, height, 20, 20);
  glPopMatrix();

  gluDeleteQuadric(quadric);
}

///////////////////////////////////////////////////////////////////////
// draw a cylinder
///////////////////////////////////////////////////////////////////////
void drawCylinder(const CYLINDER& cylinder)
{
  const VECTOR3& t = cylinder.translation();
  const Eigen::AngleAxis<GLfloat> R{ cylinder.rotation().cast<GLfloat>() };

  const REAL radius = cylinder.radius();
  const REAL height = cylinder.height();

  GLUquadricObj* quadric;
  quadric = gluNewQuadric();

  glPushMatrix();
    // apply the transform
    glTranslatef(t[0], t[1], t[2]);
    glRotatef((180.0/M_PI) * R.angle(), R.axis().x(), R.axis().y(), R.axis().z());

    // draw it along y axis instead of z axis
    glRotatef(90.0, 1.0, 0.0, 0.0);
    glTranslatef(0,0,-0.5 * height);

    // draw the end caps
    glPushMatrix();
      glTranslatef(0.0, 0.0, height);
      gluDisk(quadric, 0.0, radius, 20, 2);
    glPopMatrix();
    glPushMatrix();
      glRotatef(180.0, 0,1.0,0);
      gluDisk(quadric, 0.0, radius, 20, 2);
    glPopMatrix();

    // draw the cylinder wall
    gluCylinder(quadric, radius, radius, height, 20, 20);
  glPopMatrix();

  gluDeleteQuadric(quadric);
}

///////////////////////////////////////////////////////////////////////
// draw a sphere
///////////////////////////////////////////////////////////////////////
void drawSphere(const SPHERE& sphere)
{
  const MATRIX3& S = sphere.scale();
  const VECTOR3& t = sphere.translation();
  const Eigen::AngleAxis<GLfloat> R{ sphere.rotation().cast<GLfloat>() };

  //glColor4f(1.0, 0.0, 0.0, 0.5);
  glColor4f(0.99, 0.99, 0.99, 0.9);
  glPushMatrix();
    glTranslatef(t[0], t[1], t[2]);
    glRotatef((180.0/M_PI) * R.angle(), R.axis().x(), R.axis().y(), R.axis().z());
    glScalef(S(0,0), S(1,1), S(2,2));
    //glutSolidSphere(1.0, 20, 20);
    glutSolidSphere(1.0, 100, 100);
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////
// draw a bowl
///////////////////////////////////////////////////////////////////////
void drawBowl(const BOWL& bowl)
{
  const MATRIX3& S = bowl.scale();
  const VECTOR3& t = bowl.translation();
  const Eigen::AngleAxis<GLfloat> R{ bowl.rotation().cast<GLfloat>() };

  glCullFace(GL_FRONT);
  glColor4f(1.0, 0.0, 0.0, 0.5);
  glPushMatrix();
    // apply the transforms
    glTranslatef(t[0], t[1], t[2]);
    glRotatef((180.0/M_PI) * R.angle(), R.axis().x(), R.axis().y(), R.axis().z());
    glScalef(S(0,0), S(1,1), S(2,2));
    // Array of latitudinal triangle strips, each parallel to the equator, stacked one
    // above the other from the equator to the north pole.
    int q = 12; // subdivisions of longitudinal arc
    int p = 24; // subdivisions of latitudinal circles 
    REAL r = bowl.inRad();
    for(int j = 0; j < q; j++)
    {
      // One latitudinal triangle strip. Inner sphere
      glBegin(GL_TRIANGLE_STRIP);
        for(int i = 0; i <= p; i++)
        {
          glVertex3f( r * cos( (float)(j+1)/q * M_PI/2.0 ) * cos( 2.0 * (float)i/p * M_PI ),
                        -r * sin( (float)(j+1)/q * M_PI/2.0 ),
					    r * cos( (float)(j+1)/q * M_PI/2.0 ) * sin( 2.0 * (float)i/p * M_PI ) );
          glNormal3f( cos( (float)j/q * M_PI/2.0 ) * cos( 2.0 * (float)i/p * M_PI ),
                        -sin( (float)j/q * M_PI/2.0 ),
					    cos( (float)j/q * M_PI/2.0 ) * sin( 2.0 * (float)i/p * M_PI ) );         
          glVertex3f( r * cos( (float)j/q * M_PI/2.0 ) * cos( 2.0 * (float)i/p * M_PI ),
                        -r * sin( (float)j/q * M_PI/2.0 ),
					    r * cos( (float)j/q * M_PI/2.0 ) * sin( 2.0 * (float)i/p * M_PI ) );         
          glNormal3f( cos( (float)(j+1)/q * M_PI/2.0 ) * cos( 2.0 * (float)i/p * M_PI ),
                        -sin( (float)(j+1)/q * M_PI/2.0 ),
					    cos( (float)(j+1)/q * M_PI/2.0 ) * sin( 2.0 * (float)i/p * M_PI ) );
        }
      glEnd();
      // One latitudinal triangle strip. Outer sphere
      glBegin(GL_TRIANGLE_STRIP);
        for(int i = 0; i <= p; i++)
		    {
          glVertex3f( cos( (float)(j+1)/q * M_PI/2.0 ) * cos( 2.0 * (float)i/p * M_PI ),
                      -sin( (float)(j+1)/q * M_PI/2.0 ),
              cos( (float)(j+1)/q * M_PI/2.0 ) * sin( 2.0 * (float)i/p * M_PI ) );
          glNormal3f( cos( (float)j/q * M_PI/2.0 ) * cos( 2.0 * (float)i/p * M_PI ),
                        -sin( (float)j/q * M_PI/2.0 ),
					    cos( (float)j/q * M_PI/2.0 ) * sin( 2.0 * (float)i/p * M_PI ) );         
          glVertex3f( cos( (float)j/q * M_PI/2.0 ) * cos( 2.0 * (float)i/p * M_PI ),
                        -sin( (float)j/q * M_PI/2.0 ),
					    cos( (float)j/q * M_PI/2.0 ) * sin( 2.0 * (float)i/p * M_PI ) );         
          glNormal3f( cos( (float)(j+1)/q * M_PI/2.0 ) * cos( 2.0 * (float)i/p * M_PI ),
                        -sin( (float)(j+1)/q * M_PI/2.0 ),
					    cos( (float)(j+1)/q * M_PI/2.0 ) * sin( 2.0 * (float)i/p * M_PI ) );
		    }
      glEnd();
    }
    // Drawing latitudinal cap
    glBegin(GL_TRIANGLE_STRIP);
      for(int i = 0; i <=p; i++)
      {
        glVertex3f( r * cos( 2.0 * (float)i/p * M_PI ),
                        0, 
					    r  * sin( 2.0 * (float)i/p * M_PI ) );
        glNormal3f(0,-1,0);
        glVertex3f( cos( 2.0 * (float)i/p * M_PI ),
                       0, 
					    sin( 2.0 * (float)i/p * M_PI ) );
        glNormal3f(0,-1,0);
      }
    glEnd();

    // REFACTORING INTO JUST TRIANGLES?
  glPopMatrix();
  glCullFace(GL_BACK);
}

///////////////////////////////////////////////////////////////////////
// draw an AABB
///////////////////////////////////////////////////////////////////////
void drawAABB(const VECTOR3& minCorner, const VECTOR3& maxCorner)
{
  const VECTOR3 v000(minCorner[0], minCorner[1], minCorner[2]); 
  const VECTOR3 v100(maxCorner[0], minCorner[1], minCorner[2]); 
  const VECTOR3 v010(minCorner[0], maxCorner[1], minCorner[2]); 
  const VECTOR3 v110(maxCorner[0], maxCorner[1], minCorner[2]); 
  const VECTOR3 v001(minCorner[0], minCorner[1], maxCorner[2]); 
  const VECTOR3 v101(maxCorner[0], minCorner[1], maxCorner[2]); 
  const VECTOR3 v011(minCorner[0], maxCorner[1], maxCorner[2]); 
  const VECTOR3 v111(maxCorner[0], maxCorner[1], maxCorner[2]); 

  glColor4f(0.0, 1.0, 0.0, 0.25);
  glBegin(GL_QUADS);
    // x plus
    VECTOR3 normal = (v000 - v100).cross(v000 - v110);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v010.data());
    glVertex3dv(v110.data());
    glVertex3dv(v100.data());
    glVertex3dv(v000.data());

    // x minus
    normal = (v001 - v101).cross(v001 - v111);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v001.data());
    glVertex3dv(v101.data());
    glVertex3dv(v111.data());
    glVertex3dv(v011.data());

    // y minus
    normal = (v000 - v100).cross(v000 - v101);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v000.data());
    glVertex3dv(v100.data());
    glVertex3dv(v101.data());
    glVertex3dv(v001.data());

    // y plus
    normal = (v010 - v110).cross(v010 - v111);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v011.data());
    glVertex3dv(v111.data());
    glVertex3dv(v110.data());
    glVertex3dv(v010.data());

    // z plus
    normal = (v000 - v010).cross(v000 - v011);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v001.data());
    glVertex3dv(v011.data());
    glVertex3dv(v010.data());
    glVertex3dv(v000.data());

    // z minus
    normal = (v100 - v110).cross(v100 - v111);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v100.data());
    glVertex3dv(v110.data());
    glVertex3dv(v111.data());
    glVertex3dv(v101.data());
  glEnd();
}
void drawAABB(const AABB_NODE& node)
{
  drawAABB(node.mins, node.maxs);
}

///////////////////////////////////////////////////////////////////////
// draw a AABB tree at a specific depth
///////////////////////////////////////////////////////////////////////
void drawAABBTree(const AABB_NODE* node, const int drawDepth, const int currentDepth)
{
  if (node == NULL) return;

  if (currentDepth == drawDepth)
  {
    drawAABB(*node);
    return;
  }

  drawAABBTree(node->child[0], drawDepth, currentDepth + 1);
  drawAABBTree(node->child[1], drawDepth, currentDepth + 1);
}

///////////////////////////////////////////////////////////////////////
// draw a AABB tree at a specific depth
///////////////////////////////////////////////////////////////////////
void drawAABBTree(const AABB_TREE& tree, const int drawDepth)
{
  drawAABBTree(&(tree.root()), drawDepth, 0);
}

///////////////////////////////////////////////////////////////////////
// draw a cube
///////////////////////////////////////////////////////////////////////
void drawCube(const CUBE& cube)
{
  const MATRIX3& S = cube.scale();
  const MATRIX3& R = cube.rotation();
  const VECTOR3& t = cube.translation();
  const VECTOR3& minCorner = S * VECTOR3(-0.5, -0.5, -0.5);
  const VECTOR3& maxCorner = S * VECTOR3(0.5, 0.5, 0.5);

  VECTOR3 v000(minCorner[0], minCorner[1], minCorner[2]); 
  VECTOR3 v100(maxCorner[0], minCorner[1], minCorner[2]); 
  VECTOR3 v010(minCorner[0], maxCorner[1], minCorner[2]); 
  VECTOR3 v110(maxCorner[0], maxCorner[1], minCorner[2]); 
  VECTOR3 v001(minCorner[0], minCorner[1], maxCorner[2]); 
  VECTOR3 v101(maxCorner[0], minCorner[1], maxCorner[2]); 
  VECTOR3 v011(minCorner[0], maxCorner[1], maxCorner[2]); 
  VECTOR3 v111(maxCorner[0], maxCorner[1], maxCorner[2]); 

  v000 = R * v000 + t;
  v100 = R * v100 + t;
  v010 = R * v010 + t;
  v110 = R * v110 + t;
  v001 = R * v001 + t;
  v101 = R * v101 + t;
  v011 = R * v011 + t;
  v111 = R * v111 + t;

  glColor4f(1.0, 0.0, 0.0, 0.5);
  glBegin(GL_QUADS);
    // x plus
    VECTOR3 normal = (v000 - v100).cross(v000 - v110);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v010.data());
    glVertex3dv(v110.data());
    glVertex3dv(v100.data());
    glVertex3dv(v000.data());

    // x minus
    normal = (v001 - v101).cross(v001 - v111);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v001.data());
    glVertex3dv(v101.data());
    glVertex3dv(v111.data());
    glVertex3dv(v011.data());

    // y minus
    normal = (v000 - v100).cross(v000 - v101);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v000.data());
    glVertex3dv(v100.data());
    glVertex3dv(v101.data());
    glVertex3dv(v001.data());

    // y plus
    normal = (v010 - v110).cross(v010 - v111);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v011.data());
    glVertex3dv(v111.data());
    glVertex3dv(v110.data());
    glVertex3dv(v010.data());

    // z plus
    normal = (v000 - v010).cross(v000 - v011);
    normal.normalize();
    normal *= -1.0;
    glNormal3dv(normal.data());
    glVertex3dv(v001.data());
    glVertex3dv(v011.data());
    glVertex3dv(v010.data());
    glVertex3dv(v000.data());

    // z minus
    normal = (v100 - v110).cross(v100 - v111);
    normal.normalize();
    glNormal3dv(normal.data());
    glVertex3dv(v100.data());
    glVertex3dv(v110.data());
    glVertex3dv(v111.data());
    glVertex3dv(v101.data());
  glEnd();
}

///////////////////////////////////////////////////////////////////////
// draw coordinate axes, xyz = rgb
///////////////////////////////////////////////////////////////////////
void drawAxes()
{
  // draw coordinate axes
  glPushMatrix();
  //glTranslatef(-0.1f, -0.1f, -0.1f);
  glLineWidth(3.0f);
  glBegin(GL_LINES);
  // x axis is red
  glColor4f(10.0f, 0.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(10.0f, 0.0f, 0.0f, 0.0f);
  glVertex3f(10.0f, 0.0f, 0.0f);

  // y axis is green
  glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(0.0f, 10.0f, 0.0f, 1.0f);
  glVertex3f(0.0f, 10.0f, 0.0f);

  // z axis is blue
  glColor4f(0.0f, 0.0f, 10.0f, 1.0f);
  glVertex3f(0.0f, 0.0f, 0.0f);
  glColor4f(0.0f, 0.0f, 10.0f, 0.0f);
  glVertex3f(0.0f, 0.0f, 10.0f);
  glEnd();
  glLineWidth(1.0f);
  glPopMatrix();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawVertexFacePairs(const TET_MESH& tetMesh, const int highlighted)
{
  const vector<pair<int,int> >& vertexFacePairs = tetMesh.vertexFaceCollisions();
  const vector<VECTOR3>& vertices = tetMesh.vertices();
  const vector<VECTOR3I>& surfaceTriangles = tetMesh.surfaceTriangles();

  for (unsigned int x = 0; x < vertexFacePairs.size(); x++)
  {
    if (highlighted >= 0 && (int)x != highlighted) continue;

    // get the barycentric coordinate
    const int triangleID = vertexFacePairs[x].second;
    const int vertexID = vertexFacePairs[x].first;
    const VECTOR3I t = surfaceTriangles[triangleID];
    vector<VECTOR3> tet;
    tet.push_back(vertices[vertexID]);
    tet.push_back(vertices[t[0]]);
    tet.push_back(vertices[t[1]]);
    tet.push_back(vertices[t[2]]);

    //VECTOR3 bary = getBarycentricCoordinates(tet);
    VECTOR3 bary = getInsideBarycentricCoordinates(tet);

    if ((int)x == highlighted)
    {
      bary = getInsideBarycentricCoordinatesDebug(tet);
      cout << " vertex: " << vertexID << " triangle: " << triangleID << endl;
      //cout << " bary: " << bary.transpose() << endl;
    }

    const VECTOR3 vf = bary[0] * vertices[t[0]] + bary[1] * vertices[t[1]] + bary[2] * vertices[t[2]];

    const int v0 = vertexFacePairs[x].first;
    glBegin(GL_POINTS);
      glColor4f(10.0, 0.0, 0.0, 1.0);
      glVertex3f(vertices[v0][0], vertices[v0][1], vertices[v0][2]);
      glColor4f(0.0, 0.0, 10.0, 1.0);
      glVertex3f(vf[0], vf[1], vf[2]);
    glEnd();

    glLineWidth(5.0);
    glBegin(GL_LINES);
      glColor4f(0.0, 10.0, 0.0, 1.0);
      glVertex3f(vertices[v0][0], vertices[v0][1], vertices[v0][2]);
      glVertex3f(vf[0], vf[1], vf[2]);
    glEnd();

    //const int triangleID = vertexFacePairs[x].second;
    glColor4f(0.0, 0.0, 1.0, 0.5);
    glBegin(GL_TRIANGLES);
      glVertex3f(vertices[t[0]][0], vertices[t[0]][1], vertices[t[0]][2]);
      glVertex3f(vertices[t[1]][0], vertices[t[1]][1], vertices[t[1]][2]);
      glVertex3f(vertices[t[2]][0], vertices[t[2]][1], vertices[t[2]][2]);
    glEnd();
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawEdgeEdgePairs(const TET_MESH& tetMesh, const int pairID)
{
  const vector<pair<int,int> >& edgeEdgeCollisions = tetMesh.edgeEdgeCollisions();
  const vector<pair<VECTOR2,VECTOR2> >& edgeEdgeCoordinates = tetMesh.edgeEdgeCoordinates();
  const vector<VECTOR3>& vertices = tetMesh.vertices();
  const vector<VECTOR2I>& surfaceEdges = tetMesh.surfaceEdges();

  const vector<bool> edgeEdgeIntersections = tetMesh.edgeEdgeIntersections();

  for (unsigned int x = 0; x < edgeEdgeCollisions.size(); x++)
  {
    // if we're just drawing a specific pair, filter for that one
    if (pairID != -1 && (int)x != pairID) continue;

    vector<VECTOR3> vs(4);
    const int edge0 = edgeEdgeCollisions[x].first;
    const int edge1 = edgeEdgeCollisions[x].second;
    vs[0] = vertices[surfaceEdges[edge0][0]];
    vs[1] = vertices[surfaceEdges[edge0][1]];
    vs[2] = vertices[surfaceEdges[edge1][0]];
    vs[3] = vertices[surfaceEdges[edge1][1]];

    const VECTOR2 a = edgeEdgeCoordinates[x].first;
    const VECTOR2 b = edgeEdgeCoordinates[x].second;

    const VECTOR3 middle0 = a[0] * vs[0] + a[1] * vs[1];
    const VECTOR3 middle1 = b[0] * vs[2] + b[1] * vs[3];

    glLineWidth(5.0);

    glBegin(GL_LINES);
      glColor4f(10.0, 0.0, 0.0, 10.0);
      glVertex3f(vs[0][0], vs[0][1], vs[0][2]);
      glVertex3f(vs[1][0], vs[1][1], vs[1][2]);
      
      glColor4f(0.0, 0.0, 10.0, 10.0);
      glVertex3f(vs[2][0], vs[2][1], vs[2][2]);
      glVertex3f(vs[3][0], vs[3][1], vs[3][2]);
      
      glColor4f(0.0, 10.0, 0.0, 10.0);
      glVertex3f(middle0[0], middle0[1], middle0[2]);
      glVertex3f(middle1[0], middle1[1], middle1[2]);
    glEnd();

    // just draw the first one
    //if (x > 10) return;
  }
}

///////////////////////////////////////////////////////////////////////
void drawKinematicShape(const KINEMATIC_SHAPE& shape)
{
  using namespace std;

  const string& name = shape.name();

  if (name.compare(string("CUBE")) == 0)
  {
    drawCube((const CUBE&)shape);
    return;
  }

  if (name.compare(string("CYLINDER")) == 0)
  {
    drawCylinder((const CYLINDER&)shape);
    return;
  }
 
  if (name.compare(string("CAPSULE")) == 0)
  {
    drawCapsule((const CAPSULE&)shape);
    return;
  }

  if (name.compare(string("SPHERE")) == 0)
  {
    drawSphere((const SPHERE&)shape);
    return;
  }
  
  if (name.compare(string("BOWL")) == 0)
  {
    drawBowl((const BOWL&)shape);
    return;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawCollisions(const STRAND_MESH& mesh)
{
  const vector<pair<int, int> >& collisions = mesh.edgeEdgeCollisions();
  const vector<VECTOR2I>& edges = mesh.edgeIndices();
  const vector<VECTOR3>& vertices = mesh.vertices();
  const vector<pair<VECTOR2, VECTOR2> >& coordinates = mesh.edgeEdgeCoordinates();

  glDisable(GL_DEPTH_TEST);

  glLineWidth(10.0);
  glColor4f(0.0, 0.0, 1.0, 100.0);
  glBegin(GL_LINES);
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    const pair<int,int>& edgePair = collisions[x];
    const VECTOR2I& edge0 = edges[edgePair.first];
    const VECTOR2I& edge1 = edges[edgePair.second];

    const VECTOR2& coordinate0 = coordinates[x].first;
    const VECTOR2& coordinate1 = coordinates[x].second;

    const VECTOR3 middle0 = coordinate0[0] * vertices[edge0[0]] + coordinate0[1] * vertices[edge0[1]];
    const VECTOR3 middle1 = coordinate1[0] * vertices[edge1[0]] + coordinate1[1] * vertices[edge1[1]];
    glVertex3f(middle0[0], middle0[1], middle0[2]);
    glVertex3f(middle1[0], middle1[1], middle1[2]);
  }
  glEnd();

  glColor4f(1.0, 1.0, 0.0, 100.0);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    const pair<int,int>& edgePair = collisions[x];
    const VECTOR2I& edge0 = edges[edgePair.first];
    const VECTOR2I& edge1 = edges[edgePair.second];

    const VECTOR2& coordinate0 = coordinates[x].first;
    const VECTOR2& coordinate1 = coordinates[x].second;
    
    const VECTOR3 middle0 = coordinate0[0] * vertices[edge0[0]] + coordinate0[1] * vertices[edge0[1]];
    const VECTOR3 middle1 = coordinate1[0] * vertices[edge1[0]] + coordinate1[1] * vertices[edge1[1]];
    glVertex3f(middle0[0], middle0[1], middle0[2]);
    glVertex3f(middle1[0], middle1[1], middle1[2]);
  }
  glEnd();

  glEnable(GL_DEPTH_TEST);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void drawCollisionsOld(const STRAND_MESH& mesh)
{
  const vector<pair<int, int> >& collisions = mesh.edgeEdgeCollisionsOld();
  const vector<VECTOR2I>& edges = mesh.edgeIndices();
  const vector<VECTOR3>& vertices = mesh.verticesOld();
  const vector<pair<VECTOR2, VECTOR2> >& coordinates = mesh.edgeEdgeCoordinatesOld();

  glDisable(GL_DEPTH_TEST);

  glLineWidth(10.0);
  glColor4f(0.0, 0.0, 1.0, 100.0);
  glBegin(GL_LINES);
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    const pair<int,int>& edgePair = collisions[x];
    const VECTOR2I& edge0 = edges[edgePair.first];
    const VECTOR2I& edge1 = edges[edgePair.second];

    const VECTOR2& coordinate0 = coordinates[x].first;
    const VECTOR2& coordinate1 = coordinates[x].second;

    const VECTOR3 middle0 = coordinate0[0] * vertices[edge0[0]] + coordinate0[1] * vertices[edge0[1]];
    const VECTOR3 middle1 = coordinate1[0] * vertices[edge1[0]] + coordinate1[1] * vertices[edge1[1]];
    glVertex3f(middle0[0], middle0[1], middle0[2]);
    glVertex3f(middle1[0], middle1[1], middle1[2]);
  }
  glEnd();

  glColor4f(1.0, 1.0, 0.0, 100.0);
  glBegin(GL_POINTS);
  for (unsigned int x = 0; x < collisions.size(); x++)
  {
    const pair<int,int>& edgePair = collisions[x];
    const VECTOR2I& edge0 = edges[edgePair.first];
    const VECTOR2I& edge1 = edges[edgePair.second];

    const VECTOR2& coordinate0 = coordinates[x].first;
    const VECTOR2& coordinate1 = coordinates[x].second;

    //cout << " coordinate 0: " << coordinate0.transpose() << endl;
    //cout << " coordinate 1: " << coordinate1.transpose() << endl;
    
    const VECTOR3 middle0 = coordinate0[0] * vertices[edge0[0]] + coordinate0[1] * vertices[edge0[1]];
    const VECTOR3 middle1 = coordinate1[0] * vertices[edge1[0]] + coordinate1[1] * vertices[edge1[1]];

    //const VECTOR3 diff = middle0 - middle1;
    //cout << " collision diff: " << diff.norm() << endl;
    glVertex3f(middle0[0], middle0[1], middle0[2]);
    glVertex3f(middle1[0], middle1[1], middle1[2]);
  }
  glEnd();

  glEnable(GL_DEPTH_TEST);
}


} // HOBAK
