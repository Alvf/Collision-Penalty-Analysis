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
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "ARAP.h"
#include "util/MATRIX_UTIL.h"

namespace HOBAK {
namespace SHELL {

ARAP::ARAP(const REAL& mu, const REAL& lambda) :
  STRETCHING(mu, lambda)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string ARAP::name() const
{ 
  return std::string("ARAP"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL ARAP::psi(const MATRIX3x2& F) const
{
  MATRIX3x2 R;
  MATRIX2 S;
  polarDecomposition(F, R, S);
  return _mu * (F - R).squaredNorm();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// If nothing is provided, take the SVD and call PK1
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 ARAP::PK1(const MATRIX3x2& F) const
{
  MATRIX3x2 R;
  MATRIX2 S;
  polarDecomposition(F, R, S);

  const MATRIX2 I = MATRIX2::Identity();
  return R * (2.0 * _mu * (S - I));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// This is the eigensystem from "A finite element formulation of baraff-witkin cloth", Kim 2020
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 ARAP::hessian(const MATRIX3x2& F) const
{
  MATRIX3x2 subU;
  VECTOR2 sigma;
  MATRIX2 V;
  svd(F, subU, sigma, V);

  MATRIX3 U;
  VECTOR3 u0 = subU.col(0);
  VECTOR3 u1 = subU.col(1);

  // Panetta says this is the deformed surface normal, but
  // should be equivalent to the third direction we need.
  //
  // it doesn't need to be normalized! The two columns of U
  // are already guaranteed to be unit magnitude
  // VECTOR3 u2 = u0.cross(u1).normalized();
  VECTOR3 u2 = u0.cross(u1);
  U.col(0) = u0;
  U.col(1) = u1;
  U.col(2) = u2;

  // build out the eigenvalues
  double lambda[6];
  for (int x = 0; x < 6; x++)
    lambda[x] = 2.0;
  lambda[0] -= 4.0 / (sigma[0] + sigma[1]);
  lambda[1] -= 2.0 / sigma[0];
  lambda[2] -= 2.0 / sigma[1];

  // build out the eigenvectors
  MATRIX3x2 eigenmatrices[6];
  MATRIX3x2 middle;
  const double invSqrt2 = 1.0 / sqrt(2.0);

  middle << 0, -invSqrt2, invSqrt2, 0, 0, 0;
  eigenmatrices[0] = U * middle * V.transpose();

  middle << 0, 0, 0, 0, 1, 0;
  eigenmatrices[1] = U * middle * V.transpose();

  middle << 0, 0, 0, 0, 0, 1;
  eigenmatrices[2] = U * middle * V.transpose();

  middle << invSqrt2, 0, 0, invSqrt2, 0, 0;
  eigenmatrices[3] = U * middle * V.transpose();

  middle << invSqrt2, 0, 0, -invSqrt2, 0, 0;
  eigenmatrices[4] = U * middle * V.transpose();

  middle << 0, invSqrt2, invSqrt2, 0, 0, 0;
  eigenmatrices[5] = U * middle * V.transpose();

  MATRIX6 pPpF;
  pPpF.setZero();

  for (int x = 0; x < 6; x++)
  {
    const VECTOR6 flat = flatten(eigenmatrices[x]);
    pPpF += lambda[x] * (flat * flat.transpose());
  }

  return _mu * pPpF;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// This is the eigensystem from "A finite element formulation of baraff-witkin cloth", Kim 2020
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 ARAP::clampedHessian(const MATRIX3x2& F) const
{
  MATRIX3x2 subU;
  VECTOR2 sigma;
  MATRIX2 V;
  svd(F, subU, sigma, V);

  MATRIX3 U;
  VECTOR3 u0 = subU.col(0);
  VECTOR3 u1 = subU.col(1);

  // Panetta says this is the deformed surface normal, but
  // should be equivalent to the third direction we need.
  //
  // it doesn't need to be normalized! The two columns of U
  // are already guaranteed to be unit magnitude
  // VECTOR3 u2 = u0.cross(u1).normalized();
  VECTOR3 u2 = u0.cross(u1);
  U.col(0) = u0;
  U.col(1) = u1;
  U.col(2) = u2;

  // build out the eigenvalues
  double lambda[6];
  for (int x = 0; x < 6; x++)
    lambda[x] = 2.0;
  lambda[0] -= 4.0 / (sigma[0] + sigma[1]);
  lambda[1] -= 2.0 / sigma[0];
  lambda[2] -= 2.0 / sigma[1];

  // build out the eigenvectors
  MATRIX3x2 eigenmatrices[6];
  MATRIX3x2 middle;
  const double invSqrt2 = 1.0 / sqrt(2.0);

  middle << 0, -invSqrt2, invSqrt2, 0, 0, 0;
  eigenmatrices[0] = U * middle * V.transpose();

  middle << 0, 0, 0, 0, 1, 0;
  eigenmatrices[1] = U * middle * V.transpose();

  middle << 0, 0, 0, 0, 0, 1;
  eigenmatrices[2] = U * middle * V.transpose();

  middle << invSqrt2, 0, 0, invSqrt2, 0, 0;
  eigenmatrices[3] = U * middle * V.transpose();

  middle << invSqrt2, 0, 0, -invSqrt2, 0, 0;
  eigenmatrices[4] = U * middle * V.transpose();

  middle << 0, invSqrt2, invSqrt2, 0, 0, 0;
  eigenmatrices[5] = U * middle * V.transpose();

  MATRIX6 pPpF;
  pPpF.setZero();

  for (int x = 0; x < 6; x++)
  {
    // clamp eigenvalues
    if (lambda[x] <= 0.0) continue;

    const VECTOR6 flat = flatten(eigenmatrices[x]);
    pPpF += lambda[x] * (flat * flat.transpose());
  }

  return _mu * pPpF;
}

} // SHELL
} // HOBAK
