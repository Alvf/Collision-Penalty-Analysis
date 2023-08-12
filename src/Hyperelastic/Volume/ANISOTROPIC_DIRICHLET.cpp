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
#include "ANISOTROPIC_DIRICHLET.h"
#include "MATRIX_UTIL.h"

#define ONE_SIDED_FILTER 1

namespace HOBAK {
namespace VOLUME {

ANISOTROPIC_DIRICHLET::ANISOTROPIC_DIRICHLET(const REAL& mu, const VECTOR3& a) :
    _mu(mu), _a(a)
{
  _a.normalize();
  _A = _a * _a.transpose();
}

std::string ANISOTROPIC_DIRICHLET::name() const
{
  return "ANISOTROPIC_DIRICHLET";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL ANISOTROPIC_DIRICHLET::psi(const MATRIX3& F) const
{
#if ONE_SIDED_FILTER
  if (F.determinant() >= 0.0) return 0.0;
#endif

  return _mu * invariant5(F, _a);
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 ANISOTROPIC_DIRICHLET::PK1(const MATRIX3& F) const
{
#if ONE_SIDED_FILTER
  MATRIX3 result = 2.0 * _mu * F * _A;
  if (F.determinant() >= 0.0)
    result *= 0.0;
  return result;
#else
  return 2.0 * _mu * F * _A;
#endif
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_DIRICHLET::hessian(const MATRIX3& F) const
{
  MATRIX9 H = MATRIX9::Zero();
  
  // iterate over 3x3 blocks
  for (int y = 0; y < 3; y++)
      for (int x = 0; x < 3; x++)
      {
          const int x3 = 3 * x;
          const int y3 = 3 * y;

          for (int i = 0; i < 3; i++)
              H(x3 + i, y3 + i) = 2.0 * _mu * _A(x,y);
      }
#if ONE_SIDED_FILTER
  if (F.determinant() >= 0.0)
    H *= 0.0;
#endif

  return H;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ANISOTROPIC_DIRICHLET::clampedHessian(const MATRIX3& F) const
{
  // it can't go indefinite, so just return the Hessian
  return hessian(F);
}

bool ANISOTROPIC_DIRICHLET::energyNeedsSVD() const
{
    return false;
}

bool ANISOTROPIC_DIRICHLET::PK1NeedsSVD() const
{
    return false;
}

}
}
