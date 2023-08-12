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
#ifndef BACKWARD_EULER_VELOCITY_H
#define BACKWARD_EULER_VELOCITY_H

#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Timestepper/Volume/TIMESTEPPER.h"

namespace HOBAK {
namespace VOLUME {

////////////////////////////////////////////////////////////////////////////////////////////////////
// This is an implementation of the Baraff-Witkin-style velocity-level solver from
// "Large Steps in Cloth Simulation", SIGGRAPH 1998
//
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
////////////////////////////////////////////////////////////////////////////////////////////////////
class BACKWARD_EULER_VELOCITY : public TIMESTEPPER
{
public:
  BACKWARD_EULER_VELOCITY(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, ENERGY_12D& vfEnergy, ENERGY_12D& eeEnergy);

  // take a timestep
  virtual bool solve(const bool verbose) override;
  bool solveRayleighDamped(const bool verbose);
  bool solveEnergyDamped(const bool verbose);

private:
  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  virtual void updateConstraintTargets() override;

  // Baraff-Witkin solves for change in velocity
  VECTOR _vDelta;

  // current simulation time
  REAL _time;
 
  // current simulation step
  int _currentTimestep;
};

} // HOBAK
} // TIMESTEPPER

#endif
