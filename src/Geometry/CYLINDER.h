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
#ifndef CYLINDER_H
#define CYLINDER_H

#include "KINEMATIC_SHAPE.h"

namespace HOBAK {

// cylinder positions are defined using R * S * x + t,
// starting from a cylinder centered at (0,0,0), running up and down
// the y-axis with height 1 and radius 1
//
// Reminder: to make a new rotation matrix in Eigen, do:
// Eigen::AngleAxisd(0.1, VECTOR3::UnitX())
class CYLINDER : public KINEMATIC_SHAPE
{
public:
  CYLINDER(const VECTOR3& center, const REAL& radius, const REAL& height);
	~CYLINDER();

  const REAL radius() const { return _radius; };
  const REAL height() const { return _height; };

  virtual bool inside(const VECTOR3& point) const override;
  virtual REAL distance(const VECTOR3& point) const override;

  // for the cylinder it's slightly easier if we don't apply the scaling
  // via a matrix in the local-to-world transform
  virtual VECTOR3 localVertexToWorld(const VECTOR3& local) const override
  {
    return _rotation * local + _translation;
  };
  virtual VECTOR3 worldVertexToLocal(const VECTOR3& world) const override
  {
    return _rotation.transpose() * (world - _translation);
  };

  // remember that "inside" is negative with signed distance
  virtual REAL signedDistance(const VECTOR3& point) const override;

  // get the closest point on the cube, as well as the normal at the point
  virtual void getClosestPoint(const VECTOR3& query, 
                               VECTOR3& closestPointLocal, 
                               VECTOR3& normalLocal) const override;

protected:
  REAL _radius;
  REAL _height;
};

}

#endif
