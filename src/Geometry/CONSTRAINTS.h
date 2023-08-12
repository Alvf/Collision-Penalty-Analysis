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
#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include "SETTINGS.h"
#include "KINEMATIC_SHAPE.h"

namespace HOBAK {

// constrain a vertex to move along with the position of a kinematic shape
struct KINEMATIC_CONSTRAINT
{
  const KINEMATIC_SHAPE* shape;
  int vertexID;
  VECTOR3 localPosition;
};

// constrain a vertex to be on a plane, but can slide in the tangential directions
struct PLANE_CONSTRAINT
{
  const KINEMATIC_SHAPE* shape;
  int vertexID;
  VECTOR3 localClosestPoint;
  VECTOR3 localNormal;

  // at the end of the last timestep, was the body lifting away
  // from the surface? If so, the next time buildSurfaceConstraints()
  // is called, we should delete this constraint.
  bool isSeparating;
};

}

#endif
