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
#ifndef QUASISTATIC_H
#define QUASISTATIC_H

#include "Geometry/TET_MESH.h"
#include "Geometry/KINEMATIC_SHAPE.h"
#include "Geometry/CONSTRAINTS.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Timestepper/Volume/TIMESTEPPER.h"

namespace HOBAK {
namespace VOLUME {

class QUASISTATIC : public TIMESTEPPER
{
public:
  QUASISTATIC(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, ENERGY_12D& vfEnergy, ENERGY_12D& eeEnergy);

  // take a timestep
  //virtual bool solve(const bool verbose) override { return solve(verbose, 0.01); };
  virtual bool solve(const bool verbose) override { return solve(verbose, 1.0); };

  // take a step with a prescribed stepsize
  bool solve(const bool verbose, const REAL stepSize);

  // take a timestep, but use a backtracking line search
  bool solveWithBacktracking(const bool verbose);

private:
  // update the displacement targets the the Baraff-Witkin-style constraints
  // are trying to hit. Assumes that buildConstraintMatrix() has already been called
  virtual void updateConstraintTargets() override;

  // filter positions to incorporate Baraff-Witkin-style constraints
  virtual void applyKinematicConstraints() override;

  // find all the surface vertices that are in collision and create constraints
  virtual void findNewSurfaceConstraints(const bool verbose) override;

  int _minIterations;
  int _maxIterations;
  REAL _residualTolerance;
  
  int _seenNewtonIterations;
};

} // HOBAK
} // TIMESTEPPER

#endif
