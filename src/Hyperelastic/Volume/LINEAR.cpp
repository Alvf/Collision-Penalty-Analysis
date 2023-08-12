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
#include "LINEAR.h"
#include "MATRIX_UTIL.h"

namespace HOBAK {
namespace VOLUME {

LINEAR::LINEAR(const REAL& mu, const REAL& lambda)
: _mu(mu)
, _lambda(lambda)
{
  assert(_mu >= 0.0);
  assert(_lambda >= 0.0);
}

REAL LINEAR::psi(const MATRIX3& F) const
{
  const MATRIX3 eps = 0.5 * (F + F.transpose()) - MATRIX3::Identity();
  const REAL trEps = eps.trace();
  return _mu * eps.squaredNorm() + 0.5 * _lambda * trEps * trEps;
}

MATRIX3 LINEAR::PK1(const MATRIX3& F) const
{
  const MATRIX3 eps = 0.5 * (F + F.transpose()) - MATRIX3::Identity();
  const REAL trEps = eps.trace();
  return 2.0 * _mu * eps + _lambda * trEps * MATRIX3::Identity();
}

MATRIX9 LINEAR::hessian(const MATRIX3& F) const
{
  MATRIX9 H;
  H = _mu * MATRIX9::Identity();
  H(0,0) += _mu + _lambda;
  H(4,4) += _mu + _lambda;
  H(8,8) += _mu + _lambda;

  H(1,3) = H(3,1) = _mu;
  H(2,6) = H(6,2) = _mu;
  H(5,7) = H(7,5) = _mu;

  H(0,4) = H(0,8) = H(4,8) = _lambda;
  H(4,0) = H(8,0) = H(8,4) = _lambda;
  return H;
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
MATRIX9 LINEAR::clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return hessian(U * Sigma.asDiagonal() * V.transpose());
}

std::string LINEAR::name() const
{
  return "LINEAR";
}

bool LINEAR::energyNeedsSVD() const
{
  return false;
}

bool LINEAR::PK1NeedsSVD() const
{
  return false;
}

} // VOLUME
} // HOBAK
