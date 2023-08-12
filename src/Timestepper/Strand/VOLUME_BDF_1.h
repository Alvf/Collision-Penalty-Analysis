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
#ifndef VOLUME_BDF_1_H
#define VOLUME_BDF_1_H

#include "Timestepper/Strand/VOLUME_TIMESTEPPER.h"

namespace HOBAK {
namespace STRAND {

////////////////////////////////////////////////////////////////////////////////////////////////////
// Newton-style solver for first-order Backward Differentiation Formula (BDF-1).
//
// A multi-step version of the Baraff-Witkin solver.
////////////////////////////////////////////////////////////////////////////////////////////////////
class VOLUME_BDF_1 : public VOLUME_TIMESTEPPER
{
public:
  VOLUME_BDF_1(TET_STRAND_MESH& strandMesh, STRAND::STRETCHING& stretching, ENERGY_12D& eeEnergy);
  
  virtual bool solveDynamics(const bool verbose) override;
  virtual bool solveDynamicsWithRotation(const bool verbose, const MATRIX3& rotation, const VECTOR3& translation) override;
  
  virtual void setDt(const REAL dt) override;

private:
  REAL _alpha[6];
  VECTOR _acceleration;
  VECTOR _accelerationOld;
};

} // HOBAK
} // TIMESTEPPER

#endif
