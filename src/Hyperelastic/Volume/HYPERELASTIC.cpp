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
// May 26, 2021 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "HYPERELASTIC.h"
#include "util/MATRIX_UTIL.h"

namespace HOBAK {
namespace VOLUME {

HYPERELASTIC::~HYPERELASTIC() = default;

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, then just recombine everything into F and call that version of Psi
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL HYPERELASTIC::psi(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return psi(U * Sigma.asDiagonal() * V.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call Psi
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL HYPERELASTIC::psi(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  return psi(U, Sigma, V);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, then just recombine everything into F and call that version of PK1
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3 HYPERELASTIC::PK1(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return PK1(U * Sigma.asDiagonal() * V.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call PK1
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3 HYPERELASTIC::PK1(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  return PK1(U, Sigma, V);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, then just recombine everything into F and call that version of Hessian
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 HYPERELASTIC::hessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return hessian(U * Sigma.asDiagonal() * V.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call Hessian
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 HYPERELASTIC::hessian(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  return hessian(U, Sigma, V);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, then just recombine everything into F and call that version of 
// ClampedHessian
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 HYPERELASTIC::clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const
{
  return clampedHessian(U * Sigma.asDiagonal() * V.transpose());
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call ClampedHessian
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX9 HYPERELASTIC::clampedHessian(const MATRIX3& F) const
{
  MATRIX3 U,V;
  VECTOR3 Sigma;
  svd_rv(F, U, Sigma, V);
  return clampedHessian(U, Sigma, V);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL HYPERELASTIC::computeMu(const REAL E, const REAL nu)
{
  return E / (2.0 * (1.0 + nu));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL HYPERELASTIC::computeLambda(const REAL E, const REAL nu)
{
  return (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
}

} // VOLUME
} // HOBAK
