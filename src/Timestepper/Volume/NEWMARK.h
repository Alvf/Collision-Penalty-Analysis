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
#ifndef NEWMARK_H
#define NEWMARK_H

#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Damping/Volume/DAMPING.h"
#include "Timestepper/Volume/TIMESTEPPER.h"

namespace HOBAK {
namespace VOLUME {

////////////////////////////////////////////////////////////////////////////////////////////////////
// [TJM15] refers to "Smoothed Aggregation Multigrid for Cloth Simulation", SIGGRAPH Asia 2015
////////////////////////////////////////////////////////////////////////////////////////////////////
class NEWMARK : public TIMESTEPPER
{
public:
  NEWMARK(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, ENERGY_12D& vfEnergy, ENERGY_12D& eeEnergy);
  const REAL rayleighAlpha() const  { return _rayleighAlpha; };
  const REAL rayleighBeta() const   { return _rayleighBeta; };

  const int maxNewtonIterations() const { return _maxNewtonIterations; };
  const VECTOR velocityOld() const      { return _velocityOld; };
  VECTOR& velocityOld()      { return _velocityOld; };
  int& maxNewtonIterations() { return _maxNewtonIterations; };
  virtual void setDt(const REAL dt) override;

  // take a timestep
  virtual bool solve(const bool verbose) override;
  virtual bool solveRayleighDamped(const bool verbose);
  virtual bool solveEnergyDamped(const bool verbose);

  // output all Newton steps to OBJ files?
  // this is for just-for-rfun visualizations
  bool& outputNewton() { return _outputNewton; };

protected:
  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  virtual void updateConstraintTargets() override;

  int _minNewtonIterations;
  int _maxNewtonIterations;
  int _seenNewtonIterations;
  REAL _residualTolerance;
  
  // Newmark constants
  REAL _alpha[6];
  //REAL _accelerationAlpha[5];
  REAL _beta;
  REAL _gamma;
  
  // current simulation time
  REAL _time;
 
  // current simulation step
  int _currentTimestep;

  // variables to solve for
  VECTOR _acceleration;
  VECTOR _velocityOld;
  VECTOR _accelerationOld;

  // output all Newton steps to OBJ files?
  bool _outputNewton;
};

} // HOBAK
} // VOLUME

#endif
