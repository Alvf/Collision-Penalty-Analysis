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
#ifndef BACKWARD_EULER_POSITION_H
#define BACKWARD_EULER_POSITION_H

#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Damping/Volume/DAMPING.h"
#include "Timestepper/Volume/TIMESTEPPER.h"
#include "util/BLOCK_SPARSE_MATRIX3.h"

namespace HOBAK {
namespace VOLUME {

////////////////////////////////////////////////////////////////////////////////////////////////////
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
////////////////////////////////////////////////////////////////////////////////////////////////////
class BACKWARD_EULER_POSITION : public TIMESTEPPER
{
public:
  BACKWARD_EULER_POSITION(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, ENERGY_12D& vfEnergy, ENERGY_12D& eeEnergy);
  BACKWARD_EULER_POSITION(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, VOLUME::DAMPING& damping, ENERGY_12D& vfEnergy, ENERGY_12D& eeEnergy);
  virtual ~BACKWARD_EULER_POSITION();

  const REAL rayleighAlpha() const { return _rayleighAlpha; };
  const REAL rayleighBeta() const  { return _rayleighBeta; };

  const VECTOR velocityOld() const { return _velocityOld; };
  VECTOR& velocityOld()            { return _velocityOld; };
  virtual void setDt(const REAL dt) override { _dt = dt; };

  // take a timestep
  virtual bool solve(const bool verbose) override;
  bool solveEnergyDamped(const bool verbose);

  // solves with self collisions, and collision forces are Rayleigh damped
  bool solveRayleighDamped(const bool verbose);

private:
  // shared initialization across constructors
  void initialize();

  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  virtual void updateConstraintTargets() override;

  // current simulation time
  REAL _time;
 
  // current simulation step
  int _currentTimestep;

  // variables to solve for
  VECTOR _acceleration;
  VECTOR _velocityOld;
};

} // HOBAK
} // TIMESTEPPER

#endif
