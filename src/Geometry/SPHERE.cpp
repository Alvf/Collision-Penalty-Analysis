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
#include "SPHERE.h"
#include <iostream>

using namespace std;

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
// positions are defined using R * S * x + t
///////////////////////////////////////////////////////////////////////
SPHERE::SPHERE(const VECTOR3& center, const REAL& scale)
{
  _scale = MATRIX3::Identity() * scale; 
  _rotation = MATRIX3::Identity();
  _translation = center;
  _scaleInverse = _scale.inverse();

  _name = string("SPHERE");
}

SPHERE::~SPHERE()
{
}

///////////////////////////////////////////////////////////////////////
// is a point inside the sphere?
///////////////////////////////////////////////////////////////////////
bool SPHERE::inside(const VECTOR3& point) const 
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);
  REAL radius = transformed.norm();

  if (radius < 1.0)
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// distance to the sphere
///////////////////////////////////////////////////////////////////////
REAL SPHERE::distance(const VECTOR3& point) const 
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);
  REAL radius = transformed.norm();

  return fabs(radius - 1.0) * _scale(0,0);
}

///////////////////////////////////////////////////////////////////////
// signed distance to the sphere
// remember that "inside" is negative with signed distance
///////////////////////////////////////////////////////////////////////
REAL SPHERE::signedDistance(const VECTOR3& point) const
{
  // transform back to local coordinates
  VECTOR3 transformed = worldVertexToLocal(point);
  REAL radius = transformed.norm();

  return (radius - 1.0) * _scale(0,0);
}

//////////////////////////////////////////////////////////////////////
// get the closest point on the object, as well as the normal at 
// the point
//////////////////////////////////////////////////////////////////////
void SPHERE::getClosestPoint(const VECTOR3& query, 
                             VECTOR3& closestPointLocal, 
                             VECTOR3& normalLocal) const
{
  const VECTOR3 collisionPoint = worldVertexToLocal(query);
  closestPointLocal = collisionPoint.normalized();

  // this is the one instance where both of these are the same
  normalLocal = closestPointLocal;
}

} // HOBAK
