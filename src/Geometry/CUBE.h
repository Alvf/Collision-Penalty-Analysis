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
#ifndef CUBE_H
#define CUBE_H

#include "KINEMATIC_SHAPE.h"

namespace HOBAK {

// box positions are defined using R * S * x + t,
// starting from a cube centered at (0,0,0), with sides of length 1,
// so the max corner is at  (0.5, 0.5, 0.5)
// and the min corner is at (-0.5, -0.5, -0.5)
//
// Reminder: to make a new rotation matrix in Eigen, do:
// Eigen::AngleAxisd(0.1, VECTOR3::UnitX())
class CUBE : public KINEMATIC_SHAPE
{
public:
  CUBE(const VECTOR3& center, const REAL& scale);
	virtual ~CUBE();

  virtual bool inside(const VECTOR3& point) const override;
  virtual REAL distance(const VECTOR3& point) const override;

  // remember that "inside" is negative with signed distance
  virtual REAL signedDistance(const VECTOR3& point) const override;

  // get the closest point on the cube, as well as the normal at the point
  virtual void getClosestPoint(const VECTOR3& query, 
                               VECTOR3& closestPointLocal, 
                               VECTOR3& normalLocal) const override;
};

}

#endif
