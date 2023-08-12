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
#include "ANISOTROPIC_FUNG.h"
#include "MATRIX_UTIL.h"
#include <iostream>

namespace HOBAK {
namespace VOLUME {

using namespace std;

ANISOTROPIC_FUNG::ANISOTROPIC_FUNG(const REAL& mu, const VECTOR3& a) :
    _mu(mu), _a(a)
{
  _a.normalize();
  _A = _a * _a.transpose();
}

std::string ANISOTROPIC_FUNG::name() const
{
  return "ANISOTROPIC_FUNG";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
// Psi = mu * (exp(I5 - 1) - I5);
///////////////////////////////////////////////////////////////////////
REAL ANISOTROPIC_FUNG::psi(const MATRIX3& F) const
{
  const REAL I5 = (F * _a).squaredNorm();
  return _mu * (exp(I5 - 1.0) - I5);
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
// PK1 = mu * (exp(I5 - 1) - 1) * 2 * F * A;
///////////////////////////////////////////////////////////////////////
MATRIX3 ANISOTROPIC_FUNG::PK1(const MATRIX3& F) const
{
  const REAL I5 = (F * _a).squaredNorm();
  return 2.0 * _mu * (exp(I5 - 1.0) - 1.0) * F * _A;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
// fa = vec(F * A);
// I3 = eye(3,3);
// H5 = [A(1,1) * I3 A(1,2) * I3 A(1,3) * I3;
//       A(2,1) * I3 A(2,2) * I3 A(2,3) * I3;
//       A(3,1) * I3 A(3,2) * I3 A(3,3) * I3];
// H = 2 * mu * ((exp(I5 - 1) - 1) * H5 + 2 * exp(I5 - 1) * fa * fa');
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_FUNG::hessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  const VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);
  const REAL I5minusOne = I5 - 1.0;
  MATRIX9 pPpF = (exp(I5minusOne) - 1.0) * H5 + 2.0 * exp(I5minusOne) * (fa * fa.transpose());
  pPpF *= 2.0 * _mu;

  return pPpF;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_FUNG::clampedHessian(const MATRIX3& F) const
{
  const MATRIX9 H5 = kronIdentity(_A);
  const MATRIX3 FA = F * _A;
  const VECTOR9 fa = flatten(FA);

  const REAL I5 = invariant5(F, _a);
  const REAL lambda1 =           2.0 * _mu * (exp(I5 - 1.0) - 1.0);
  const REAL lambda0 = lambda1 + 4.0 * _mu * I5 * exp(I5 - 1.0);

  // try the case where it's all positive
  if (lambda0 > 0.0 && lambda1 > 0.0)
  {
    const REAL I5minusOne = I5 - 1.0;
    MATRIX9 pPpF = (exp(I5minusOne) - 1.0) * H5 + 2.0 * exp(I5minusOne) * (fa * fa.transpose());
    pPpF *= 2.0 * _mu;
    return pPpF;
  }

  // if this one's negative, you're out of luck, just return all zeros
  if (lambda0 < 0.0) return MATRIX9::Zero();

  // the only possibility left is that a rank-one subspace is positive definite
  const VECTOR9 faNormed = fa.normalized();
  return lambda0 * (faNormed * faNormed.transpose());
}

bool ANISOTROPIC_FUNG::energyNeedsSVD() const
{
  return true;
}

bool ANISOTROPIC_FUNG::PK1NeedsSVD() const
{
  return true;
}

} // VOLUME
} // HOBAK
