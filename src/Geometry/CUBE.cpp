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
#include "CUBE.h"
#include <iostream>

using namespace std;

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
// box positions are defined using R * S * x + t
///////////////////////////////////////////////////////////////////////
CUBE::CUBE(const VECTOR3& center, const REAL& scale) 
{
  _scale = MATRIX3::Identity() * scale;
  _rotation = MATRIX3::Identity();
  _translation = center;
  _scaleInverse = _scale.inverse();

  _name = string("CUBE");
}

CUBE::~CUBE()
{
}

///////////////////////////////////////////////////////////////////////
// is a point inside the box?
///////////////////////////////////////////////////////////////////////
bool CUBE::inside(const VECTOR3& point) const 
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);

  if (transformed[0] <= 0.5 && transformed[0] >= -0.5 &&
      transformed[1] <= 0.5 && transformed[1] >= -0.5 &&
      transformed[2] <= 0.5 && transformed[2] >= -0.5)
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// distance to the box
///////////////////////////////////////////////////////////////////////
REAL CUBE::distance(const VECTOR3& point) const 
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);

  if (inside(point))
  {
    REAL xMin = std::min((0.5 - transformed[0]), (transformed[0] - (-0.5)));
    REAL yMin = std::min((0.5 - transformed[1]), (transformed[1] - (-0.5)));
    REAL zMin = std::min((0.5 - transformed[2]), (transformed[2] - (-0.5)));

    return fabs(std::min(std::min(xMin, yMin), zMin));
  }
 
  // handle edges and points
  VECTOR3 diff = VECTOR3::Zero();
  if (transformed[2] > 0.5) 
    diff[2] = transformed[2] - 0.5;
  else if (transformed[2] < -0.5)
    diff[2] = -0.5 - transformed[2];
  
  if (transformed[1] > 0.5)
    diff[1] = transformed[1] - 0.5;
  else if (transformed[1] < -0.5)
    diff[1] = -0.5 - transformed[1];
  
  if (transformed[0] > 0.5)
    diff[0] = transformed[0] - 0.5;
  else if (transformed[0] < -0.5)
    diff[0] = -0.5 - transformed[0];
      
  return diff.norm() * _scale(0,0);
}

///////////////////////////////////////////////////////////////////////
// signed distance to the box
///////////////////////////////////////////////////////////////////////
REAL CUBE::signedDistance(const VECTOR3& point) const
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);

  if (inside(point))
  {
    REAL xMin = std::min((0.5 - transformed[0]), (transformed[0] - (-0.5)));
    REAL yMin = std::min((0.5 - transformed[1]), (transformed[1] - (-0.5)));
    REAL zMin = std::min((0.5 - transformed[2]), (transformed[2] - (-0.5)));

    return -fabs(std::min(std::min(xMin, yMin), zMin)) * _scale(0,0);
  }
 
  // handle edges and points
  VECTOR3 diff = VECTOR3::Zero();
  if (transformed[2] > 0.5) 
    diff[2] = transformed[2] - 0.5;
  else if (transformed[2] < -0.5)
    diff[2] = -0.5 - transformed[2];
  
  if (transformed[1] > 0.5)
    diff[1] = transformed[1] - 0.5;
  else if (transformed[1] < -0.5)
    diff[1] = -0.5 - transformed[1];
  
  if (transformed[0] > 0.5)
    diff[0] = transformed[0] - 0.5;
  else if (transformed[0] < -0.5)
    diff[0] = -0.5 - transformed[0];
      
  return diff.norm() * _scale(0,0);
}

//////////////////////////////////////////////////////////////////////
// get the closest point on the cube, as well as the normal at the point
//////////////////////////////////////////////////////////////////////
void CUBE::getClosestPoint(const VECTOR3& query, VECTOR3& closestPoint, VECTOR3& normal) const
{
  const VECTOR3 collisionPoint = worldVertexToLocal(query);

  VECTOR diffs(6);
  diffs[0] = 0.5 + collisionPoint[0];
  diffs[1] = 0.5 - collisionPoint[0];

  diffs[2] = 0.5 + collisionPoint[1];
  diffs[3] = 0.5 - collisionPoint[1];

  diffs[4] = 0.5 + collisionPoint[2];
  diffs[5] = 0.5 - collisionPoint[2];

  int minIndex = 0;
  REAL minFound = diffs[0];

  for (int x = 1; x < 6; x++)
  {
    if (diffs[x] < minFound)
    {
      minFound = diffs[x];
      minIndex = x;
    }
  }
  closestPoint = collisionPoint;
  closestPoint[0] =  -0.5;
  normal = VECTOR3(-1, 0, 0);

  switch (minIndex)
  {
    case 1:
      closestPoint = collisionPoint;
      closestPoint[0] = 0.5;
      normal = VECTOR3(1,0, 0);
      break;
    case 2:
      closestPoint = collisionPoint;
      closestPoint[1] = -0.5;
      normal = VECTOR3(0, -1,0);
      break;
    case 3:
      closestPoint = collisionPoint;
      closestPoint[1] = 0.5;
      normal = VECTOR3(0, 1,0);
      break;
    case 4:
      closestPoint = collisionPoint;
      closestPoint[2] = -0.5;
      normal = VECTOR3(0, 0, -1);
      break;
    case 5:
      closestPoint = collisionPoint;
      closestPoint[2] = 0.5;
      normal = VECTOR3(0, 0, 1);
      break;
  }
}

} // HOBAK
