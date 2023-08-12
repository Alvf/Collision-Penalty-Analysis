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
#include "DIRICHLET.h"
#include "MATRIX_UTIL.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

DIRICHLET::DIRICHLET(const REAL& mu) :
    _mu(mu)
{
}

std::string DIRICHLET::name() const
{
  return "DIRICHLET";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL DIRICHLET::psi(const MATRIX3& F) const
{
  return _mu * F.squaredNorm();
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 DIRICHLET::PK1(const MATRIX3& F) const
{
  return 2.0 * _mu * F;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 DIRICHLET::hessian(const MATRIX3& F) const
{
  MATRIX9 H;
  H.setIdentity();
  H *= 2.0 * _mu;

  return H;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 DIRICHLET::clampedHessian(const MATRIX3& F) const
{
  // it can't go indefinite, so just return the Hessian
  return hessian(F);
}

bool DIRICHLET::energyNeedsSVD() const
{
    return false;
}

bool DIRICHLET::PK1NeedsSVD() const
{
    return false;
}

}
}
