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
#include "ARAP_DAMPING.h"
#include "MATRIX_UTIL.h"
#include "../../Hyperelastic/Volume/HYPERELASTIC.h"

namespace HOBAK {
namespace VOLUME {

ARAP_DAMPING::ARAP_DAMPING(const REAL& mu)
{
  _mu = mu;
}

std::string ARAP_DAMPING::name() const
{
  return "ARAP_DAMPING";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL ARAP_DAMPING::psi(const MATRIX3& F, const MATRIX3& Fdot) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  MATRIX3 Rdot = rotationDot(U, Sigma, V, Fdot);
  return _mu * ddot((Fdot - Rdot), Fdot);
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 ARAP_DAMPING::PK1(const MATRIX3& F, const MATRIX3& Fdot) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  MATRIX3 Rdot = rotationDot(U, Sigma, V, Fdot);
  return 2.0 * _mu * (Fdot - Rdot);
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. Fdot
///////////////////////////////////////////////////////////////////////
MATRIX9 ARAP_DAMPING::hessian(const MATRIX3& F, const MATRIX3& Fdot) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  MATRIX9 pPpF = 2.0 * _mu * (MATRIX9::Identity() - rotationGradient(U,Sigma,V));

  return pPpF;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. Fdot
///////////////////////////////////////////////////////////////////////
MATRIX9 ARAP_DAMPING::clampedHessian(const MATRIX3& F, const MATRIX3& Fdot) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);

  MATRIX9 Q;
  buildTwistAndFlipEigenvectors(U,V,Q);
  buildScalingEigenvectors(U,V,Q);

  const REAL& s0 = Sigma[0];
  const REAL& s1 = Sigma[1];
  const REAL& s2 = Sigma[2];

  // init everything to 2 * mu
  VECTOR9 lambda = VECTOR9::Constant(2.0 * _mu);

  // the first few get modified
  lambda[0] = _mu * (2.0  - 4.0 / (s1 + s2));
  lambda[1] = _mu * (2.0  - 4.0 / (s0 + s2));
  lambda[2] = _mu * (2.0  - 4.0 / (s0 + s1));

  for (int x = 0; x < 9; x++)
    lambda[x] = (lambda[x] >= 0.0) ? lambda[x] : 0.0;

  return Q * lambda.asDiagonal() * Q.transpose();
}

}
}
