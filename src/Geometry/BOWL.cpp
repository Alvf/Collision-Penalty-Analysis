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
#include "BOWL.h"
#include "util/FILE_IO.h"
#include <iostream>
#include <cmath>

using namespace std;

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
// positions are defined using R * S * x + t
///////////////////////////////////////////////////////////////////////
BOWL::BOWL(const VECTOR3& center, const REAL& inRadius, const REAL& scale)
{
  _scale = MATRIX3::Identity() * scale; 
  _rotation = MATRIX3::Identity();
  _translation = center;
  _scaleInverse = _scale.inverse();

  _name = string("BOWL");

  _inRadius = inRadius;
}

BOWL::~BOWL()
{
}

///////////////////////////////////////////////////////////////////////
// is a point inside the BOWL?
///////////////////////////////////////////////////////////////////////
bool BOWL::inside(const VECTOR3& point) const 
{
  const VECTOR3 collisionPoint = worldVertexToLocal(point);

  if (collisionPoint[1] <= 0.0){
    REAL rad = collisionPoint.norm();
    return (rad <= 1 && rad >= _inRadius); 
  }
  return false;
}

///////////////////////////////////////////////////////////////////////
// distance to the BOWL
///////////////////////////////////////////////////////////////////////
REAL BOWL::distance(const VECTOR3& point) const 
{
  const VECTOR3 collisionPoint = worldVertexToLocal(point);

  if (collisionPoint[1] <= 0.0){
    REAL rad = collisionPoint.norm();
    if (rad >= 1){
      return rad - 1;
    }
    if (rad >= _inRadius){
      return fmin(1-rad, fmin(rad-_inRadius, -collisionPoint[1]));
    }
    else{
      return _inRadius - rad;
    }
  }
  else{
    const VECTOR3 cVec(collisionPoint[0],0,collisionPoint[2]); 
    REAL cRad2 = cVec.squaredNorm(); 
    if (cRad2 >= 1){
      return (collisionPoint - cVec.normalized()).norm();
      }
    else if (cRad2 <= _inRadius*_inRadius){
      return (collisionPoint - cVec.normalized()*_inRadius).norm();
    }
    else{
      return collisionPoint[1];
    }
  }
}

///////////////////////////////////////////////////////////////////////
// signed distance to the BOWL
// remember that "inside" is negative with signed distance
///////////////////////////////////////////////////////////////////////
REAL BOWL::signedDistance(const VECTOR3& point) const
{
  const VECTOR3 collisionPoint = worldVertexToLocal(point);

  if (collisionPoint[1] <= 0.0){
    REAL rad = collisionPoint.norm();
    if (rad >= 1){
      return rad - 1;
    }
    if (rad >= _inRadius){
      return fmax(rad-1, fmax(_inRadius - rad, collisionPoint[1]));
    }
    else{
      return _inRadius - rad;
    }
  }
  else{
    const VECTOR3 cVec(collisionPoint[0],0,collisionPoint[2]); 
    REAL cRad2 = cVec.squaredNorm(); 
    if (cRad2 >= 1){
      return (collisionPoint - cVec.normalized()).norm();
      }
    else if (cRad2 <= _inRadius*_inRadius){
      return (collisionPoint - cVec.normalized()*_inRadius).norm();
    }
    else{
      return collisionPoint[1];
    }
  }
}

//////////////////////////////////////////////////////////////////////
// get the closest point on the object, as well as the normal at 
// the point
//////////////////////////////////////////////////////////////////////
void BOWL::getClosestPoint(const VECTOR3& query, 
                             VECTOR3& closestPointLocal, 
                             VECTOR3& normalLocal) const
{
  const VECTOR3 collisionPoint = worldVertexToLocal(query);

  if (collisionPoint[1] <= 0.0){
    REAL rad = collisionPoint.norm();
    if (rad >= 1){
      closestPointLocal = collisionPoint.normalized();
      normalLocal = closestPointLocal;
    }
    else if (rad >= _inRadius){
      REAL din, dout, dtop;
      din = rad - _inRadius;
      dout = 1 - rad;
      dtop = -collisionPoint[1];
      REAL closest = fmin(din, fmin(dout, dtop));
      if (closest == din){
        normalLocal = -collisionPoint.normalized();
        closestPointLocal = -normalLocal*_inRadius;
      }
      else if (closest == dout){
        closestPointLocal = collisionPoint.normalized();
        normalLocal = closestPointLocal;
      }
      else if (closest == dtop){
        closestPointLocal = VECTOR3(collisionPoint[0],0,collisionPoint[2]);
        normalLocal = VECTOR3(0,1,0);
      }
    }
    else{
        normalLocal = -collisionPoint.normalized();
        closestPointLocal = -normalLocal*_inRadius;
    }
  }
  else{
    const VECTOR3 cVec(collisionPoint[0],0,collisionPoint[2]); 
    REAL cRad2 = cVec.squaredNorm(); 
    if (cRad2 >= 1){
      closestPointLocal = cVec.normalized();
      normalLocal = (collisionPoint - closestPointLocal).normalized();
    }
    else if (cRad2 <= _inRadius*_inRadius){
      closestPointLocal = cVec.normalized()*_inRadius;
      normalLocal = (collisionPoint - closestPointLocal).normalized();
    }
    else{
      closestPointLocal = cVec;
      normalLocal = VECTOR3(0,1,0);
    }
  }

}

void BOWL::writeToObj(string fileName, int p, int q) const
{
  const REAL latiStep = 2.0 * M_PI / (REAL) p;
  const REAL longStep = M_PI / (REAL) q / 2.0;
  const REAL outerR = _scale(0,0);
  const REAL innerR = _inRadius * _scale(0,0);
  const VECTOR3& center = _translation;
  vector<VECTOR3> vertices;
  vector<VECTOR3I> faces; // zero indexed

  // outer
  // first layer
  const VECTOR3 bottomOut = center + VECTOR3(0.0, -outerR, 0.0);
  int bottomIdx = 0;
  vertices.push_back(bottomOut);
  REAL layerY = - outerR * cos(longStep);
  REAL layerR = outerR * sin(longStep);
  REAL layerX = 0.0, layerZ = 0.0;
  for(int i = 0; i < p; i++){
    layerX = layerR * cos(latiStep * i);
    layerZ = -layerR * sin(latiStep * i);
    vertices.push_back(center + VECTOR3(layerX, layerY, layerZ));
    faces.push_back(VECTOR3I(bottomIdx, bottomIdx + 1 + (i + 1)%p, bottomIdx + 1 + i));
  }
  // upper layers
  for(int j = 1; j < q; j ++){
    layerY = - outerR * cos(longStep * (j + 1));
    layerR = outerR * sin(longStep * (j + 1));
    for(int i = 0; i < p; i++){
      layerX = layerR * cos(latiStep * i);
      layerZ = -layerR * sin(latiStep * i);
      vertices.push_back(center + VECTOR3(layerX, layerY, layerZ));
      faces.push_back(VECTOR3I(bottomIdx + 1 + (j - 1) * p + i, 
                               bottomIdx + 1 + j * p + (i + 1)%p, 
                               bottomIdx + 1 + j * p + i));
      faces.push_back(VECTOR3I(bottomIdx + 1 + (j - 1) * p + i,
                               bottomIdx + 1 + (j - 1) * p + (i + 1)%p,
                               bottomIdx + 1 + j * p + (i + 1)%p));
    }
  }

  // inner
  // first layer
  const VECTOR3 bottomIn = center + VECTOR3(0.0, -innerR, 0.0);
  bottomIdx = 1 + p * q;
  vertices.push_back(bottomIn);
  layerY = - innerR * cos(longStep);
  layerR = innerR * sin(longStep);
  for(int i = 0; i < p; i++){
    layerX = layerR * cos(latiStep * i);
    layerZ = -layerR * sin(latiStep * i);
    vertices.push_back(center + VECTOR3(layerX, layerY, layerZ));
    faces.push_back(VECTOR3I(bottomIdx, bottomIdx + 1 + i, bottomIdx + 1 + (i + 1)%p));
  }
  // upper layers
  for(int j = 1; j < q; j ++){
    layerY = - innerR * cos(longStep * (j + 1));
    layerR = innerR * sin(longStep * (j + 1));
    for(int i = 0; i < p; i++){
      layerX = layerR * cos(latiStep * i);
      layerZ = -layerR * sin(latiStep * i);
      vertices.push_back(center + VECTOR3(layerX, layerY, layerZ));
      faces.push_back(VECTOR3I(bottomIdx + 1 + (j - 1) * p + i, 
                               bottomIdx + 1 + j * p + i,
                               bottomIdx + 1 + j * p + (i + 1)%p));
      faces.push_back(VECTOR3I(bottomIdx + 1 + (j - 1) * p + i,
                               bottomIdx + 1 + j * p + (i + 1)%p,
                               bottomIdx + 1 + (j - 1) * p + (i + 1)%p));
    }
  }

  // brim
  const int outer0 = 1 + (q - 1) * p;
  const int inner0 = 2 + (2 * q - 1) * p;
  for(int i = 0; i < p; i++) {
    faces.push_back(VECTOR3I(outer0 + i, outer0 + (i + 1)%p, inner0 + i));
    faces.push_back(VECTOR3I(inner0 + i, outer0 + (i + 1)%p, inner0 + (i + 1)%p));
  }

  cout<< "Writing bowl to " << fileName << "..."<<endl;

  writeOBJFile(fileName.c_str(), vertices, faces);
  cout<< "done."<<endl;
}

} // HOBAK
