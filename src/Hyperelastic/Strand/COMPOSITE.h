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
///////////////////////////////////////////////////////////////////////////////////////////////////
// This is a file in the HOBAK library
//
// November 21, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef STRAND_COMPOSITE_H
#define STRAND_COMPOSITE_H

#include "SETTINGS.h"
#include "MATRIX_UTIL.h"
#include <array>

#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/STRETCHING.h"
#include "Hyperelastic/Volume/TET_STRAND_TWIST.h"

namespace HOBAK {
namespace STRAND {

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////////////
// This is a container that holds all the hyperelastic energies for a single tet-strand
///////////////////////////////////////////////////////////////////////////////////////////////////
class COMPOSITE
{
public:
  COMPOSITE();
  virtual ~COMPOSITE();

  // Computes the strain energy density
  virtual REAL psi() const;

  // The name of the material
  virtual std::string name() const;

  // accessors
  const vector<VOLUME::HYPERELASTIC*>& volumeEnergies() const { return _volumeEnergies; };
  vector<VOLUME::HYPERELASTIC*>& volumeEnergies()     { return _volumeEnergies; };
  std::array<ISOTROPIC_BENDING*,2>& bendingEnergies() { return _bendingEnergies; };
  std::array<STRETCHING*,3>& stretchingEnergies()     { return _stretchingEnergies; };

  // get all the volumetric terms (except twist)
  MATRIX3 volumePK1(const MATRIX3& F);
  MATRIX9 volumeClampedHessian(const MATRIX3& F);
  
  // get the terms for Alvin's twist energy
  MATRIX3 twistPK1(const MATRIX3& F);
  MATRIX9 twistClampedHessian(const MATRIX3& F);
  MATRIX12 twistClampedForceGradient(const vector<VECTOR3>& vertices);
  REAL twistEnergy(const MATRIX3& F);

  // get all the edge-based terms
  std::array<VECTOR3,3> edgePK1(const std::array<VECTOR3,3>& fs);
  std::array<MATRIX3,3> edgeClampedHessian(const std::array<VECTOR3,3>& fs);
  std::array<MATRIX3,3> edgeHessian(const std::array<VECTOR3,3>& fs);
  std::array<REAL,3> edgeEnergies(const std::array<VECTOR3,3>& fs);
  
  // get all the bend-based terms
  std::array<MATRIX3x2,2> bendPK1(const std::array<MATRIX3x2,2>& Es);
  std::array<MATRIX6,2> bendClampedHessian(const std::array<MATRIX3x2,2>& Es);
  std::array<REAL,2> bendEnergies(const std::array<MATRIX3x2,2>& Es);

  std::array<MATRIX3x2,2> bendPK1(const std::array<MATRIX3x2,2>& Es, const bool inverted);
  std::array<MATRIX6,2> bendClampedHessian(const std::array<MATRIX3x2,2>& Es, const bool inverted);

  // print out debugging info?
  bool& verbose() { return _verbose; };

protected:
  vector<VOLUME::HYPERELASTIC*> _volumeEnergies;

  // entry 0 is the bending between vertices 0,1,2
  // entry 1 is the bending between vertices 1,2,3
  std::array<ISOTROPIC_BENDING*,2> _bendingEnergies;

  // entry 0 is stretching between vertices 0 and 1
  // entry 1 is stretching between vertices 1 and 2
  // entry 2 is stretching between vertices 2 and 3
  std::array<STRETCHING*,3> _stretchingEnergies;

  bool _verbose;
};

}
}

#endif
