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
#include "ANISOTROPIC_STVK.h"
#include "MATRIX_UTIL.h"

// DEBUG
#include <iostream>

namespace HOBAK {
namespace VOLUME {

using namespace std;

ANISOTROPIC_STVK::ANISOTROPIC_STVK(const REAL& mu, const VECTOR3& a) :
    _mu(mu), _a(a)
{
  _a.normalize();
  _A = _a * _a.transpose();
}

std::string ANISOTROPIC_STVK::name() const
{
  return "ANISOTROPIC_STVK";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL ANISOTROPIC_STVK::psi(const MATRIX3& F) const
{
  const REAL I5minusOne = (F * _a).squaredNorm() - 1.0;
  return 0.5 * _mu * I5minusOne * I5minusOne;
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 ANISOTROPIC_STVK::PK1(const MATRIX3& F) const
{
  const REAL I5minusOne = (F * _a).squaredNorm() - 1.0;
  return (2.0 * _mu * I5minusOne) * (F * _A);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_STVK::hessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);
  MATRIX9 pPpF = (I5 - 1.0) * H5 + 2.0 * fa * fa.transpose();
  pPpF *= 2.0 * _mu;

  return pPpF;
}


///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_STVK::clampedHessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  const VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);

  const REAL lambda1 =           2.0 * _mu * (I5 - 1.0);
  const REAL lambda0 = lambda1 + 4.0 * _mu * I5;

  // try the case where it's all positive
  if (lambda0 > 0.0 && lambda1 > 0.0)
  {
    MATRIX9 pPpF = (I5 - 1.0) * H5 + 2.0 * fa * fa.transpose();
    pPpF *= 2.0 * _mu;
    return pPpF;
  }

  // if this one's negative, you're out of luck, just return all zeros
  if (lambda0 < 0.0) return MATRIX9::Zero();

  // the only possibility left is that a rank-one subspace is positive definite
  const VECTOR9 faNormed = fa.normalized();
  return lambda0 * (faNormed * faNormed.transpose());
}

bool ANISOTROPIC_STVK::energyNeedsSVD() const
{
  return true;
}

bool ANISOTROPIC_STVK::PK1NeedsSVD() const
{
  return true;
}

} // VOLUME
} // HOBAK
