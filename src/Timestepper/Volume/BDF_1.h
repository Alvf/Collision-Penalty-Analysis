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
#ifndef BDF_1_H
#define BDF_1_H

#include "Timestepper/Volume/NEWMARK.h"

namespace HOBAK {
namespace VOLUME {

////////////////////////////////////////////////////////////////////////////////////////////////////
// Newton-style solver for first-order Backward Differentiation Formula (BDF-1).
//
// A multi-step version of the Baraff-Witkin solver. It can be written as a special case
// of the Newmark equations, which is why it inherits from that parent class
////////////////////////////////////////////////////////////////////////////////////////////////////
class BDF_1 : public NEWMARK
{
public:
  BDF_1(TET_MESH& tetMesh, VOLUME::HYPERELASTIC& hyperelastic, ENERGY_12D& vfEnergy, ENERGY_12D& eeEnergy);
  virtual void setDt(const REAL dt) override;

  // take a timestep
  virtual bool solveRayleighDamped(const bool verbose) override;
  virtual bool solveEnergyDamped(const bool verbose) override;
};

} // HOBAK
} // TIMESTEPPER

#endif
