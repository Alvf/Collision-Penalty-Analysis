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

#include "QUADRATIC_STRETCHING.h"
#include "util/MATRIX_UTIL.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace STRAND {

QUADRATIC_STRETCHING::QUADRATIC_STRETCHING(const REAL& mu) :
  STRETCHING(mu)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string QUADRATIC_STRETCHING::name() const
{ 
  return std::string("QUADRATIC_STRETCHING"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL QUADRATIC_STRETCHING::psi(const VECTOR3& f) const
{
  const REAL I = f.norm();
  return _mu * 0.5 * (I - 1.0) * (I - 1.0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 QUADRATIC_STRETCHING::PK1(const VECTOR3& f) const
{
  const REAL I = f.norm();
  const REAL Iinv = 1.0 / I;

  return (_mu * (I - 1.0) * Iinv ) * f;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// position-based gradient
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR6 QUADRATIC_STRETCHING::spatialGradient(const std::vector<VECTOR3>& p, const REAL& dmInv) const
{
  const VECTOR3 f = (p[1] - p[0]) * dmInv;
  const REAL I = f.norm();
  const REAL Iinv = 1.0 / I;
  const VECTOR3 P = dmInv * (_mu * (I - 1.0) * Iinv) * f;

  VECTOR6 result;
  result.block<3,1>(0,0) = -P;
  result.block<3,1>(3,0) = P;
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3 QUADRATIC_STRETCHING::hessian(const VECTOR3& f) const
{
  const REAL I = f.norm();

  const REAL dPsi  = _mu * (I - 1.0);
  const REAL dPsi2 = _mu;

  const REAL Iinv = 1.0 / I;
  return Iinv * Iinv * (dPsi2 - Iinv * dPsi) * (f * f.transpose()) +
         Iinv * dPsi * MATRIX3::Identity();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_STRETCHING::spatialHessian(const std::vector<VECTOR3>& p, const REAL& dmInv) const
{
  const VECTOR3 f = (p[1] - p[0]) * dmInv;

  const REAL I = f.norm();
  const REAL Iinv = 1.0 / I;

  const REAL dPsi  = _mu * (I - 1.0);
  const REAL dPsi2 = _mu;

  const VECTOR3 q0 = f * Iinv;
  const REAL lambda0 = dPsi2;
  const REAL lambda1 = Iinv * dPsi;

  MATRIX2 A;
  A << 1, -1, -1, 1; 
  const MATRIX6 kron = kronIdentity(A);

  VECTOR6 qx;
  const VECTOR3 scaled = q0 * dmInv;
  qx.block<3,1>(0,0) = -scaled;
  qx.block<3,1>(3,0) = scaled;

  // this could be optimized by only computing the upper triangle of the outer-product
  // and mirroring it, and then populating the kron term by hand, since it's mostly zeros
  return (lambda0 - lambda1) * qx * qx.transpose() + (lambda1 * dmInv * dmInv) * kron;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_STRETCHING::spatialClampedHessian(const std::vector<VECTOR3>& p, const REAL& dmInv) const
{
  const VECTOR3 f = (p[1] - p[0]) * dmInv;

  const REAL I = f.norm();
  const REAL Iinv = 1.0 / I;

  const REAL dPsi  = _mu * (I - 1.0);
  const REAL dPsi2 = _mu;

  const VECTOR3 q0 = f * Iinv;
  const REAL lambda0 = dPsi2;
  const REAL lambda1 = (dPsi > 0.0) ? Iinv * dPsi : 0.0;

  MATRIX2 A;
  A << 1, -1, -1, 1; 
  const MATRIX6 kron = kronIdentity(A);

  VECTOR6 qx;
  const VECTOR3 scaled = q0 * dmInv;
  qx.block<3,1>(0,0) = -scaled;
  qx.block<3,1>(3,0) = scaled;

  // this could be optimized by only computing the upper triangle of the outer-product
  // and mirroring it, and then populating the kron term by hand, since it's mostly zeros
  return (lambda0 - lambda1) * qx * qx.transpose() + (lambda1 * dmInv * dmInv) * kron;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 QUADRATIC_STRETCHING::spatialGaussNewtonHessian(const std::vector<VECTOR3>& p, const REAL& dmInv) const
{
  const VECTOR3 f = (p[1] - p[0]) * dmInv;

  const REAL I = f.norm();
  const REAL Iinv = 1.0 / I;

  const REAL dPsi  = _mu * (I - 1.0);
  const REAL dPsi2 = _mu;

  const VECTOR3 q0 = f * Iinv;
  const REAL lambda0 = dPsi2;
  const REAL lambda1 = (dPsi > 0.0) ? Iinv * dPsi : 0.0;

  MATRIX2 A;
  A << 1, -1, -1, 1; 
  const MATRIX6 kron = kronIdentity(A);

  VECTOR6 qx;
  const VECTOR3 scaled = q0 * dmInv;
  qx.block<3,1>(0,0) = -scaled;
  qx.block<3,1>(3,0) = scaled;

  // just do outer product term
  return (lambda0 - lambda1) * qx * qx.transpose();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3 QUADRATIC_STRETCHING::clampedHessian(const VECTOR3& f) const
{
  const REAL I = f.norm();
  const REAL Iinv = 1.0 / I;

  const REAL dPsi  = _mu * (I - 1.0);
  const REAL dPsi2 = _mu;

  const VECTOR3 q0 = f * Iinv;
  const REAL lambda0 = dPsi2;
  const REAL lambda1 = (dPsi > 0.0) ? Iinv * dPsi : 0.0;

  return (lambda0 - lambda1) * (q0 * q0.transpose()) + lambda1 * MATRIX3::Identity();
}

} // STRAND 
} // HOBAK
