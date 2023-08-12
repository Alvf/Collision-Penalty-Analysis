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
#include "ANISOTROPIC_ARAP.h"
#include "MATRIX_UTIL.h"

#include <iostream>

using namespace std;

#define USING_I4 1

namespace HOBAK {
namespace VOLUME {

using namespace std;

ANISOTROPIC_ARAP::ANISOTROPIC_ARAP(const REAL& mu, const VECTOR3& a) :
    _mu(mu), _a(a)
{
  _a.normalize();
  _A = _a * _a.transpose();
}

std::string ANISOTROPIC_ARAP::name() const
{
  return "ANISOTROPIC_ARAP";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL ANISOTROPIC_ARAP::psi(const MATRIX3& F) const
{
  const REAL I5 = invariant5(F, _a);
#if USING_I4
  const REAL I4 = invariant4(F, _a);
#else
  const REAL I4 = invariant3(F);
#endif
  const REAL SI4 = (I4 > 0.0) ? 1.0 : -1.0;

  const REAL diff = sqrt(I5) - SI4;
  return 0.5 * _mu * diff * diff;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 ANISOTROPIC_ARAP::PK1(const MATRIX3& F) const
{
  const REAL I5 = invariant5(F, _a);
#if USING_I4
  const REAL I4 = invariant4(F, _a);
#else
  const REAL I4 = invariant3(F);
#endif
  const REAL SI4 = (I4 > 0.0) ? 1.0 : -1.0;

  // if we're not near the singularity, just return the generic
  // formula
  const REAL diracEps = (sizeof(REAL) == sizeof(double)) ? 1e-8 : 1e-4;
  if (fabs(I4) >= diracEps)
    return _mu * (1.0 - SI4 / sqrt(I5)) * (F * _A);

  return -1.0 * _mu * (1.0 - 1.0 / 2.0) * (F * _A);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_ARAP::hessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);
  const REAL sqrtI5inv = 1.0 / sqrt(I5);
#if USING_I4
  const REAL I4 = invariant4(F, _a);
#else
  const REAL I4 = invariant3(F);
#endif
  const REAL SI4 = (I4 > 0.0) ? 1.0 : -1.0;

  const MATRIX9 pPpF = (1.0  - SI4 * sqrtI5inv) * H5 + 
                       (SI4 * sqrtI5inv * sqrtI5inv * sqrtI5inv) * (fa  * fa.transpose());

  return _mu * pPpF;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
//
// this assumes that _mu is positive, the subtractive stiffness
// described in the paper isn't here yet
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_ARAP::clampedHessian(const MATRIX3& F) const
{
  const REAL I5 = invariant5(F, _a);
#if USING_I4
  const REAL I4 = invariant4(F, _a);
#else
  const REAL I4 = invariant3(F);
#endif
  const REAL SI4 = (I4 > 0.0) ? 1.0 : -1.0;

  // if it is not under compression, the Hessian is guaranteed to
  // be positive definite
  //if (1.0 - SI4 / sqrt(I5) >= 0.0)
  if (I5 >= SI4)
    return hessian(F);
 
  const MATRIX3 FA = F * _A;
  const VECTOR9 fa = flatten(FA);
  return _mu * (1.0 / I5) * (fa * fa.transpose());
}

bool ANISOTROPIC_ARAP::energyNeedsSVD() const
{
  return true;
}

bool ANISOTROPIC_ARAP::PK1NeedsSVD() const
{
  return true;
}

} // VOLUME
} // HOBAK
