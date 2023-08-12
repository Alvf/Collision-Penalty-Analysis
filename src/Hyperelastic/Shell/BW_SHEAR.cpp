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
// It is the shearing energy from the original Baraff-Witkin paper,
// "Large Steps in Cloth Simulation" SIGGRAPH 1998
//
// But, it's the FEM-style formulation and eigenanalysis from the Kim paper
// "A Finite Element Formulation of Baraff-Witkin Cloth" SCA 2020
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "BW_SHEAR.h"
#include "util/MATRIX_UTIL.h"

namespace HOBAK {
namespace SHELL {

BW_SHEAR::BW_SHEAR(const REAL& mu, const REAL& lambda) :
  STRETCHING(mu, lambda)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string BW_SHEAR::name() const
{ 
  return std::string("BW_SHEAR"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL BW_SHEAR::psi(const MATRIX3x2& F) const
{
  const VECTOR2 u(1.0, 0.0);
  const VECTOR2 v(0.0, 1.0);
  const REAL I6 = (F * u).transpose() * (F * v);

  return _mu * 0.5 * I6 * I6;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call PK1
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 BW_SHEAR::PK1(const MATRIX3x2& F) const
{
  const VECTOR2 u(1.0, 0.0);
  const VECTOR2 v(0.0, 1.0);
  const REAL I6 = (F * u).transpose() * (F * v);
  MATRIX3x2 P = F * (u * v.transpose() + v * u.transpose());

  // half from mu and leading 2 on P cancel
  P *= _mu * I6;

  return P;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// This is the eigensystem from "A finite element formulation of baraff-witkin cloth", Kim 2020
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 BW_SHEAR::hessian(const MATRIX3x2& F) const
{
  const VECTOR2 u(1.0, 0.0);
  const VECTOR2 v(0.0, 1.0);
  const REAL I6 = (F * u).transpose() * (F * v);

  const MATRIX3x2 product = F * (u * v.transpose() + v * u.transpose());
  const VECTOR6 gradient = flatten(product);

  MATRIX6 pPpF;
  pPpF.setZero();

  // lower-left block-identity
  pPpF(3, 0) = I6;
  pPpF(4, 1) = I6;
  pPpF(5, 2) = I6;

  // upper-right block-identity
  pPpF(0, 3) = I6;
  pPpF(1, 4) = I6;
  pPpF(2, 5) = I6;

  pPpF += gradient * gradient.transpose();

  // half from mu and leading 2 on P cancel
  pPpF *= _mu;

  return pPpF;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// This is the eigensystem from "A finite element formulation of baraff-witkin cloth", Kim 2020
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 BW_SHEAR::clampedHessian(const MATRIX3x2& F) const
{
  const VECTOR2 u(1.0, 0.0);
  const VECTOR2 v(0.0, 1.0);
  const REAL I6 = (F * u).transpose() * (F * v);
  const REAL signI6 = (I6 >= 0) ? 1.0 : -1.0;

  const MATRIX3x2 product = F * (u * v.transpose() + v * u.transpose());
  const VECTOR6 gradient = flatten(product);

  MATRIX6 H;
  H.setZero();
  H(3,0) = H(4,1) = H(5,2) = H(0,3) = H(1,4) = H(2,5) = 1.0;

  // get the novel eigenvalue
  const REAL I2 = F.squaredNorm();
  const REAL lambda0 = 0.5 * (I2 + sqrt(I2 * I2 + 12.0 * I6 * I6));

  // get the novel eigenvector
  // the H matrix multiply is a column swap; could be optimized further
  const VECTOR6 q0 = (I6 * H * gradient + lambda0 * gradient).normalized();

  MATRIX6 T = MATRIX6::Identity();
  T = 0.5 * (T + signI6 * H);

  const VECTOR6 Tq = T * q0;
  const REAL normTq = Tq.squaredNorm();

  MATRIX6 pPpF = fabs(I6) * (T - (Tq * Tq.transpose()) / normTq) + lambda0 * (q0 * q0.transpose());

  // half from mu and leading 2 on Hessian cancel
  pPpF *= _mu;

  return pPpF;
}

} // SHELL
} // HOBAK
