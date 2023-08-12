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
#include "ARAP.h"
#include "MATRIX_UTIL.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

ARAP::ARAP(const REAL& mu, const REAL& lambda ) :
    _lambda(lambda), _mu(mu)
{
}

std::string ARAP::name() const
{
  return "ARAP";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy, using the SVD
///////////////////////////////////////////////////////////////////////
REAL ARAP::psi(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  const MATRIX3 F = U * Sigma.asDiagonal() * V.transpose();
  return _mu * (F - U * V.transpose()).squaredNorm();
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 ARAP::PK1(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  MATRIX3 R = U * V.transpose();
  MATRIX3 S = V * Sigma.asDiagonal() * V.transpose();

  const MATRIX3 I = MATRIX3::Identity();
  return R * (2.0 * _mu * (S - I));
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
///////////////////////////////////////////////////////////////////////
MATRIX9 ARAP::hessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  const REAL& s0 = Sigma[0];
  const REAL& s1 = Sigma[1];
  const REAL& s2 = Sigma[2];

  // create the pseudo-twist vectors
  MATRIX3 twist0, twist1, twist2; 
  twist0 <<  0, 0, 0,
             0, 0, -1,
             0, 1, 0;
  twist1 <<  0, 0, 1,
             0, 0, 0,
            -1, 0, 0;
  twist2 <<  0, 1, 0,
            -1, 0, 0,
             0, 0, 0;

  // compute the eigenvectors
  const REAL front = 1.0 / sqrt(2.0);
  const MATRIX3 Q0 = front * U * twist0 * V.transpose();
  const MATRIX3 Q1 = front * U * twist1 * V.transpose();
  const MATRIX3 Q2 = front * U * twist2 * V.transpose();

  // flatten them out to vectors
  const VECTOR9 q0 = flatten(Q0);
  const VECTOR9 q1 = flatten(Q1);
  const VECTOR9 q2 = flatten(Q2);

  // compute the eigenvectors
  const REAL lambda0 = _mu * 2.0 / (s1 + s2);
  const REAL lambda1 = _mu * 2.0 / (s0 + s2);
  const REAL lambda2 = _mu * 2.0 / (s0 + s1);

  // populate the trivial eigenvectors
  MATRIX9 pPpF;
  pPpF.setIdentity();
  pPpF *= _mu;

  // put in the non-trivial ones
  pPpF -= lambda0 * (q0 * q0.transpose());
  pPpF -= lambda1 * (q1 * q1.transpose());
  pPpF -= lambda2 * (q2 * q2.transpose());
  pPpF *= 2.0;
  return pPpF;
}

///////////////////////////////////////////////////////////////////////
// derivative of PK1 w.r.t. F
//
// this is from Section 5.1 in
// "Analytic Eigensystems for Isotropic Distortion Energies"
///////////////////////////////////////////////////////////////////////
MATRIX9 ARAP::clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
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

bool ARAP::energyNeedsSVD() const
{
    return true;
}

bool ARAP::PK1NeedsSVD() const
{
    return true;
}

} // VOLUME
} // HOBAK
