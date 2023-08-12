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

#include "STVK.h"
#include "util/MATRIX_UTIL.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace SHELL {

STVK::STVK(const REAL& mu, const REAL& lambda) :
  STRETCHING(mu, lambda)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string STVK::name() const
{ 
  return std::string("STVK"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL STVK::psi(const MATRIX3x2& F) const
{
  //const MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  //return _mu * E.squaredNorm() + _lambda * 0.5 * pow(E.trace(), 2.0); 

  // invariant-based version
  const REAL I2 = invariant2(F);
  const REAL I3 = invariant3(F);
  return _mu * 0.25 * (I2 * I2 - 2.0 * I3 * I3 - 2.0 * I2 + 2.0) +
         (_lambda / 8.0) * (I2 - 2.0) * (I2 - 2.0);

  /*
  cout << " diff: " << test - E.squaredNorm() << endl;
  cout << " test: " << test << endl;
  cout << " ground: " << E.squaredNorm() << endl;
  cout << " ratio: " << test / E.squaredNorm() << endl;
  const REAL ground = pow(E.trace(), 2.0);
  cout << " diff: " << test2 - ground << endl;
  cout << " test: " << test2 << endl;
  cout << " ground: " << ground << endl;
  cout << endl;
  */

}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 STVK::PK1(const MATRIX3x2& F) const
{
  // build the second PK
  const MATRIX2 E = 0.5 * (F.transpose() * F - MATRIX2::Identity());
  MATRIX2 S = _lambda * E.trace() * MATRIX2::Identity() + 2.0 * _mu * E;

  // convert to PK1
  return F * S;
}

static MATRIX3x2 partialFij(const int i, const int j)
{
  MATRIX3x2 A;
  A.setZero();
  A(i,j) = 1.0;
  return A;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 STVK::hessian(const MATRIX3x2& F) const
{
  MATRIX6 stretch;
  stretch.setZero();
  MATRIX6 volume;
  volume.setZero();
  
  int index = 0;
  for (int y = 0; y < 2; y++)
    for (int x = 0; x < 3; x++, index++)
    {
      MATRIX3x2 Fij = partialFij(x,y);
      MATRIX3x2 A = 4.0 * (Fij * F.transpose() * F + F * Fij.transpose() * F +
                           F * F.transpose() * Fij - Fij);
      stretch.col(index) = flatten(A);

      MATRIX3x2 B = (F.squaredNorm() - 2.0) * Fij + 2.0 * F(x,y) * F;
      volume.col(index) = flatten(B);
    }
  return 0.25 * (_mu * stretch + 2.0 * _lambda * volume);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 STVK::clampedHessian(const MATRIX3x2& F) const
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

  const REAL I2 = invariant2(F);
  const REAL I3 = invariant3(F);

  // build out the eigenvalues
  VECTOR6 lambda;
  //lambda[0] = 2.0 * _mu * sigma[0] * sigma[0] - _mu * sigma[1] * sigma[1] + _mu * (I2 - 1);
  //lambda[1] = 2.0 * _mu * sigma[1] * sigma[1] - _mu * sigma[0] * sigma[0] + _mu * (I2 - 1);
  lambda[2] = _mu * (I2 - 1) - _mu * I3 + _lambda * (0.5 * I2 - 1.0);
  lambda[3] = _mu * (I2 - 1) + _mu * I3 + _lambda * (0.5 * I2 - 1.0);
  lambda[4] = _mu * (I2 - 1) - _mu * I3 * (sigma[1] / sigma[0]) + _lambda * (0.5 * I2 - 1.0);
  lambda[5] = _mu * (I2 - 1) - _mu * I3 * (sigma[0] / sigma[1]) + _lambda * (0.5 * I2 - 1.0);

  // build out the system matrix
  MATRIX2 A; 
  A(0,0) = 4.0 * sigma[0] * sigma[0] * (_lambda / 4.0 + _mu / 2.0) - _mu * sigma[1] * sigma[1] + _mu * (I2 - 1) + _lambda * (0.5 * I2 - 1.0);
  A(1,1) = 4.0 * sigma[1] * sigma[1] * (_lambda / 4.0 + _mu / 2.0) - _mu * sigma[0] * sigma[0] + _mu * (I2 - 1) + _lambda * (0.5 * I2 - 1.0);
  A(0,1) = -2.0 * _mu * I3 + 4.0 * I3 * (_lambda / 4.0 + _mu  / 2.0);
  A(1,0) = A(0,1);

  MATRIX2 Q;
  VECTOR2 Lambda;
  eigensystem(A, Q, Lambda);

  lambda[0] = Lambda[0];
  lambda[1] = Lambda[1];

  // build out the eigenvectors
  MATRIX3x2 eigenmatrices[6];
  MATRIX3x2 middle;
  const double invSqrt2 = 1.0 / sqrt(2.0);

  VECTOR2 q = Q.col(0);
  middle << q[0], 0, 0, q[1], 0, 0;
  eigenmatrices[0] = U * middle * V.transpose();

  q = Q.col(1);
  middle << q[0], 0, 0, q[1], 0, 0;
  eigenmatrices[1] = U * middle * V.transpose();

  middle << 0, -invSqrt2, invSqrt2, 0, 0, 0;
  eigenmatrices[2] = U * middle * V.transpose();

  middle << 0, invSqrt2, invSqrt2, 0, 0, 0;
  eigenmatrices[3] = U * middle * V.transpose();

  middle << 0, 0, 0, 0, 1, 0;
  eigenmatrices[4] = U * middle * V.transpose();

  middle << 0, 0, 0, 0, 0, 1;
  eigenmatrices[5] = U * middle * V.transpose();

  MATRIX6 pPpF;
  pPpF.setZero();
  for (int x = 0; x < 6; x++)
  {
    if (lambda[x] <= 0.0) continue;

    const VECTOR6 flat = flatten(eigenmatrices[x]);
    pPpF += lambda[x] * (flat * flat.transpose());
  }

  return pPpF;
}

} // SHELL
} // HOBAK
