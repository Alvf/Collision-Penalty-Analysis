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
///////////////////////////////////////////////////////////////////////
VECTOR3 planeNormal(const vector<VECTOR3>& plane)
{
  const VECTOR3 edge1 = plane[1] - plane[0];
  const VECTOR3 edge2 = plane[2] - plane[0];
  return edge1.cross(edge2).normalized();
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

} // HOBAK
