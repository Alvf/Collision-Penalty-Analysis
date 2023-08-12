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
// May 26, 2022 Theodore Kim, kim@cs.yale.edu
//
// It is the stretching energy from the original Baraff-Witkin paper,
// "Large Steps in Cloth Simulation" SIGGRAPH 1998
//
// But, it's the FEM-style formulation and eigenanalysis from the Kim paper
// "A Finite Element Formulation of Baraff-Witkin Cloth" SCA 2020
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "BW_STRETCH.h"
#include "util/MATRIX_UTIL.h"

namespace HOBAK {
namespace SHELL {

BW_STRETCH::BW_STRETCH(const REAL& mu, const REAL& lambda) :
  STRETCHING(mu, lambda)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string BW_STRETCH::name() const
{ 
  return std::string("BW_STRETCH"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL BW_STRETCH::psi(const MATRIX3x2& F) const
{
  const VECTOR3 wu = F.col(0);
  const VECTOR3 wv = F.col(1);

  VECTOR2 C;
  //C << wu.norm() - _bu, wv.norm() - _bv;
  C << wu.norm() - 1.0, wv.norm() - 1.0;
  return _mu * 0.5 * (C.dot(C));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call PK1
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 BW_STRETCH::PK1(const MATRIX3x2& F) const
{
  const VECTOR2 u(1.0, 0.0);
  const VECTOR2 v(0.0, 1.0);
  const REAL I5u = (F * u).transpose() * (F * u);
  const REAL I5v = (F * v).transpose() * (F * v);
  const REAL invSqrtI5u = 1.0 / sqrt(I5u);
  const REAL invSqrtI5v = 1.0 / sqrt(I5v);
  const MATRIX2 uOuter = u * u.transpose();
  const MATRIX2 vOuter = v * v.transpose();

  MATRIX3x2 P = F - (invSqrtI5u * F * uOuter + invSqrtI5v * F * vOuter);

  // half from mu and leading 2 on P cancel
  P *= _mu;

  return P;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// This is the eigensystem from "A finite element formulation of baraff-witkin cloth", Kim 2020
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 BW_STRETCH::hessian(const MATRIX3x2& F) const
{
  MATRIX6 pPpF;
  pPpF.setZero();
  const VECTOR2 u(1.0, 0.0);
  const VECTOR2 v(0.0, 1.0);
  const REAL I5u = (F * u).transpose() * (F * u);
  const REAL I5v = (F * v).transpose() * (F * v);
  const REAL invSqrtI5u = 1.0 / sqrt(I5u);
  const REAL invSqrtI5v = 1.0 / sqrt(I5v);

  // set the block diagonals,
  pPpF(0, 0) = pPpF(1, 1) = pPpF(2, 2) = (1.0 - invSqrtI5u);
  pPpF(3, 3) = pPpF(4, 4) = pPpF(5, 5) = (1.0 - invSqrtI5v);

  // upper block diagonal,
  const VECTOR3 fu = F.col(0).normalized();
  const REAL uCoeff = invSqrtI5u;
  pPpF.block<3,3>(0,0) += uCoeff * (fu * fu.transpose());

  // lower block diagonal
  const VECTOR3 fv = F.col(1).normalized();
  const REAL vCoeff = invSqrtI5v;
  pPpF.block<3,3>(3,3) += vCoeff * (fv * fv.transpose());

  // the leading 2 is absorbed by the mu / 2 coefficient
  pPpF *= _mu;

  return pPpF;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// This is the eigensystem from "A finite element formulation of baraff-witkin cloth", Kim 2020
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 BW_STRETCH::clampedHessian(const MATRIX3x2& F) const
{
  MATRIX6 pPpF;
  pPpF.setZero();
  const VECTOR2 u(1.0, 0.0);
  const VECTOR2 v(0.0, 1.0);
  const REAL I5u = (F * u).transpose() * (F * u);
  const REAL I5v = (F * v).transpose() * (F * v);
  const REAL invSqrtI5u = 1.0 / sqrt(I5u);
  const REAL invSqrtI5v = 1.0 / sqrt(I5v);

  // set the block diagonals,
  // build the rank-three subspace with all-(1 / invSqrtI5) eigenvalues
  pPpF(0, 0) = pPpF(1, 1) = pPpF(2, 2) = std::max((1.0 - invSqrtI5u), 0.0);
  pPpF(3, 3) = pPpF(4, 4) = pPpF(5, 5) = std::max((1.0 - invSqrtI5v), 0.0);

  // modify the upper block diagonal,
  // bump the single rank-one eigenvector back to just 1, unless it
  // was clamped, then just set it directly to 1
  const VECTOR3 fu = F.col(0).normalized();
  const REAL uCoeff = (1.0 - invSqrtI5u >= 0.0) ? invSqrtI5u : 1.0;
  pPpF.block<3,3>(0,0) += uCoeff * (fu * fu.transpose());

  // modify the lower block diagonal similarly
  const VECTOR3 fv = F.col(1).normalized();
  const REAL vCoeff = (1.0 - invSqrtI5v >= 0.0) ? invSqrtI5v : 1.0;
  pPpF.block<3,3>(3,3) += vCoeff * (fv * fv.transpose());

  // the leading 2 is absorbed by the mu / 2 coefficient
  pPpF *= _mu;

  return pPpF;
}

} // SHELL
} // HOBAK
