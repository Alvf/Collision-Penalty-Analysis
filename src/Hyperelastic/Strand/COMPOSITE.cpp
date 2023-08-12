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
// July 11, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "COMPOSITE.h"
#include "Hyperelastic/Volume/ANISOTROPIC_ARAP.h"
#include "Hyperelastic/Volume/TET_STRAND_TWIST.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace STRAND {

COMPOSITE::COMPOSITE()
{
  for (unsigned int x = 0; x < 2; x++)
    _bendingEnergies[x] = NULL;
  
  for (unsigned int x = 0; x < 3; x++)
    _stretchingEnergies[x] = NULL;

  _verbose = false;
}

COMPOSITE::~COMPOSITE()
{
  for (unsigned int x = 0; x < _volumeEnergies.size(); x++)
    delete _volumeEnergies[x];
  for (unsigned int x = 0; x < 2; x++)
    delete _bendingEnergies[x];
  for (unsigned int x = 0; x < 3; x++)
    delete _stretchingEnergies[x];
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string COMPOSITE::name() const
{ 
  return std::string("COMPOSITE"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL COMPOSITE::psi() const
{
  return 0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL COMPOSITE::twistEnergy(const MATRIX3& F)
{
  REAL energy = 0;
  for (unsigned int x = 0; x < _volumeEnergies.size(); x++)
  {
    // this is somewhat inefficient, but unlikely to be a bottleneck
    if (_volumeEnergies[x]->name().compare("TET_STRAND_TWIST") != 0)
      continue;

    energy += _volumeEnergies[x]->psi(F);
  }

  return energy;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3 COMPOSITE::twistPK1(const MATRIX3& F)
{
  MATRIX3 PK1;
  PK1.setZero();
  for (unsigned int x = 0; x < _volumeEnergies.size(); x++)
  {
    // this is somewhat inefficient, but unlikely to be a bottleneck
    if (_volumeEnergies[x]->name().compare("TET_STRAND_TWIST") != 0)
      continue;

    PK1 += _volumeEnergies[x]->PK1(F);
  }

  return PK1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3 COMPOSITE::volumePK1(const MATRIX3& F)
{
  MATRIX3 PK1;
  PK1.setZero();
  for (unsigned int x = 0; x < _volumeEnergies.size(); x++)
  {
    // this is somewhat inefficient, but unlikely to be a bottleneck
    if (_volumeEnergies[x]->name().compare("TET_STRAND_TWIST") == 0)
      continue;

    PK1 += _volumeEnergies[x]->PK1(F);
  }

  /*
  if (_verbose)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    cout << " F: " << endl << F << endl;
    cout << " det(F): " << F.determinant() << endl;
    cout << " volume terms: " << _volumeEnergies.size() << endl;
    cout << " Psi: " << _volumeEnergies[0]->psi(F) << endl;

    const VOLUME::ANISOTROPIC_ARAP* arap = (VOLUME::ANISOTROPIC_ARAP*)_volumeEnergies[0];
    const VECTOR3 a = arap->anisotropyDirection();
    cout << " a: " << a.transpose() << endl;
  
    const REAL I5 = invariant5(F, a);
    const REAL I4 = invariant4(F, a);
    const REAL SI4 = (I4 > 0.0) ? 1.0 : -1.0;
    const REAL diff = sqrt(I5) - SI4;
    cout << " sqrt(I5):   " << sqrt(I5) << endl;
    cout << " I4:         " << I4 << endl;
    cout << " SI4:        " << SI4 << endl;
    cout << " diff:       " << diff << endl;
  }
  */

  return PK1;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 COMPOSITE::twistClampedHessian(const MATRIX3& F)
{
  MATRIX9 H;
  H.setZero();
  for (unsigned int x = 0; x < _volumeEnergies.size(); x++)
  {
    // this is somewhat inefficient, but unlikely to be a bottleneck
    if (_volumeEnergies[x]->name().compare("TET_STRAND_TWIST") != 0)
      continue;
   
    H += _volumeEnergies[x]->clampedHessian(F);
  }

  return H;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX12 COMPOSITE::twistClampedForceGradient(const vector<VECTOR3>& vertices)
{
  MATRIX12 H;
  H.setZero();
  for (unsigned int x = 0; x < _volumeEnergies.size(); x++)
  {
    // this is somewhat inefficient, but unlikely to be a bottleneck
    if (_volumeEnergies[x]->name().compare("TET_STRAND_TWIST") != 0)
      continue;

    using namespace HOBAK::VOLUME;

    TET_STRAND_TWIST* twist = (TET_STRAND_TWIST*)_volumeEnergies[x];
    H += twist->clampedForceGradient(vertices);
  }

  return H;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 COMPOSITE::volumeClampedHessian(const MATRIX3& F)
{
  MATRIX9 H;
  H.setZero();
  for (unsigned int x = 0; x < _volumeEnergies.size(); x++)
  {
    // this is somewhat inefficient, but unlikely to be a bottleneck
    if (_volumeEnergies[x]->name().compare("TET_STRAND_TWIST") == 0)
      continue;
    
    H += _volumeEnergies[x]->clampedHessian(F);
  }

  return H;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the edge-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<REAL,3> COMPOSITE::edgeEnergies(const std::array<VECTOR3,3>& fs)
{
  std::array<REAL,3> energies;
  for (unsigned int x = 0; x < 3; x++)
  {
    energies[x] = 0;

    if (_stretchingEnergies[x] != NULL)
      energies[x] = _stretchingEnergies[x]->psi(fs[x]);
  }

  return energies;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the edge-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<VECTOR3,3> COMPOSITE::edgePK1(const std::array<VECTOR3,3>& fs)
{
  std::array<VECTOR3,3> PK1s;
  for (unsigned int x = 0; x < 3; x++)
  {
    PK1s[x].setZero();

    if (_stretchingEnergies[x] != NULL)
      PK1s[x] = _stretchingEnergies[x]->PK1(fs[x]);
  }

  return PK1s;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the bend-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<REAL,2> COMPOSITE::bendEnergies(const std::array<MATRIX3x2,2>& Es)
{
  std::array<REAL,2> energy;
  for (unsigned int x = 0; x < 2; x++)
  {
    energy[x] = 0;  

    if (_bendingEnergies[x] != NULL)
      energy[x] = _bendingEnergies[x]->psi(Es[x]);
  }

  return energy;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the bend-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<MATRIX3x2,2> COMPOSITE::bendPK1(const std::array<MATRIX3x2,2>& Es, const bool inverted)
{
  std::array<MATRIX3x2,2> PK1s;
  for (unsigned int x = 0; x < 2; x++)
  {
    PK1s[x].setZero();

    if (_bendingEnergies[x] != NULL)
      PK1s[x] = _bendingEnergies[x]->PK1(Es[x], inverted);
  }

  if (_verbose)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    cout << " Energy: " << _bendingEnergies[0]->psi(Es[0], inverted) << endl;
  }

  return PK1s;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the bend-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<MATRIX3x2,2> COMPOSITE::bendPK1(const std::array<MATRIX3x2,2>& Es)
{
  std::array<MATRIX3x2,2> PK1s;
  for (unsigned int x = 0; x < 2; x++)
  {
    PK1s[x].setZero();

    if (_bendingEnergies[x] != NULL)
      PK1s[x] = _bendingEnergies[x]->PK1(Es[x]);
  }
  
  if (_verbose)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    cout << " Energy: " << _bendingEnergies[0]->psi(Es[0]) << endl;
    cout << " PK1 norm: " << PK1s[0].norm() << endl;

    const REAL theta = ISOTROPIC_BENDING::angle(Es[0]);
    const REAL theta0 = _bendingEnergies[0]->theta0();
    const REAL diff = (theta - theta0);

    cout << " theta0: " << theta0 << endl;
    cout << " theta:  " << theta << endl;
    cout << " (theta - theta0)^2: " << diff * diff << endl;
    cout << " Theta (degrees): " << ISOTROPIC_BENDING::angle(Es[0]) * (360.0 / (2.0 * M_PI)) << endl;
    cout << " E norm: " << endl << Es[0].norm() << endl;
    cout << " E: " << endl << Es[0] << endl;
  }

  return PK1s;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the edge-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<MATRIX6,2> COMPOSITE::bendClampedHessian(const std::array<MATRIX3x2,2>& Es, const bool inverted)
{
  std::array<MATRIX6,2> Hs;
  for (unsigned int x = 0; x < 2; x++)
  {
    Hs[x].setZero();

    if (_bendingEnergies[x] != NULL)
      Hs[x] = _bendingEnergies[x]->clampedHessian(Es[x], inverted);
  }

  return Hs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the edge-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<MATRIX6,2> COMPOSITE::bendClampedHessian(const std::array<MATRIX3x2,2>& Es)
{
  std::array<MATRIX6,2> Hs;
  for (unsigned int x = 0; x < 2; x++)
  {
    Hs[x].setZero();

    if (_bendingEnergies[x] != NULL)
      Hs[x] = _bendingEnergies[x]->clampedHessian(Es[x]);
  }

  return Hs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the edge-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<MATRIX3,3> COMPOSITE::edgeClampedHessian(const std::array<VECTOR3,3>& fs)
{
  std::array<MATRIX3,3> Hs;
  for (unsigned int x = 0; x < 3; x++)
  {
    Hs[x].setZero();

    if (_stretchingEnergies[x] != NULL)
      Hs[x] = _stretchingEnergies[x]->clampedHessian(fs[x]);
  }

  return Hs;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// get all the edge-based terms
///////////////////////////////////////////////////////////////////////////////////////////////////
std::array<MATRIX3,3> COMPOSITE::edgeHessian(const std::array<VECTOR3,3>& fs)
{
  std::array<MATRIX3,3> Hs;
  for (unsigned int x = 0; x < 3; x++)
  {
    Hs[x].setZero();

    if (_stretchingEnergies[x] != NULL)
      Hs[x] = _stretchingEnergies[x]->hessian(fs[x]);
  }

  return Hs;
}

} // STRAND 
} // HOBAK
