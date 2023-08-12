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

#include "THETA_FASTER.h"
#include "util/MATRIX_UTIL.h"
#include "TIMER.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace SHELL {

THETA_FASTER::THETA_FASTER(const REAL& mu) :
  THETA(mu)
{
  initializeCrossHessians();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string THETA_FASTER::name() const
{ 
  return std::string("THETA_FASTER"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// populate _crossHessian0 and _crossHessian1
///////////////////////////////////////////////////////////////////////////////////////////////////
void THETA_FASTER::initializeCrossHessians()
{
  for (int y = 0; y < 12; y++)
    for (int x = 0; x < 12; x++)
    {
      _crossHessian0[x][y].setZero();
      _crossHessian1[x][y].setZero();
    }

  _crossHessian0[0][4] = _crossHessian0[4][0] =
  _crossHessian0[1][6] = _crossHessian0[6][1] =
  _crossHessian0[3][7] = _crossHessian0[7][3] = VECTOR3(0,0,-1);
  
  _crossHessian0[0][5] = _crossHessian0[5][0] =
  _crossHessian0[2][6] = _crossHessian0[6][2] =
  _crossHessian0[3][8] = _crossHessian0[8][3] = VECTOR3(0,1,0);
  
  _crossHessian0[0][7] = _crossHessian0[7][0] =
  _crossHessian0[1][3] = _crossHessian0[3][1] =
  _crossHessian0[4][6] = _crossHessian0[6][4] = VECTOR3(0,0,1);
 
  _crossHessian0[0][8] = _crossHessian0[8][0] =
  _crossHessian0[2][3] = _crossHessian0[3][2] =
  _crossHessian0[5][6] = _crossHessian0[6][5] = VECTOR3(0,-1,0);

  _crossHessian0[1][5] = _crossHessian0[5][1] =
  _crossHessian0[2][7] = _crossHessian0[7][2] =
  _crossHessian0[4][8] = _crossHessian0[8][4] = VECTOR3(-1,0,0);
  
  _crossHessian0[1][8] = _crossHessian0[8][1] =
  _crossHessian0[2][4] = _crossHessian0[4][2] =
  _crossHessian0[5][7] = _crossHessian0[7][5] = VECTOR3(1,0,0);

  _crossHessian1[3][7] = _crossHessian1[7][3] =
  _crossHessian1[4][9] = _crossHessian1[9][4] =
  _crossHessian1[6][10] = _crossHessian1[10][6] = VECTOR3(0,0,1);

  _crossHessian1[3][8] = _crossHessian1[8][3] =
  _crossHessian1[6][11] = _crossHessian1[11][6] =
  _crossHessian1[5][9] = _crossHessian1[9][5] = VECTOR3(0,-1,0);

  _crossHessian1[7][9] = _crossHessian1[9][7] =
  _crossHessian1[4][6] = _crossHessian1[6][4] =
  _crossHessian1[3][10] = _crossHessian1[10][3] = VECTOR3(0,0,-1);

  _crossHessian1[8][9] = _crossHessian1[9][8] =
  _crossHessian1[5][6] = _crossHessian1[6][5] =
  _crossHessian1[3][11] = _crossHessian1[11][3] = VECTOR3(0,1,0);

  _crossHessian1[4][8] = _crossHessian1[8][4] =
  _crossHessian1[7][11] = _crossHessian1[11][7] =
  _crossHessian1[5][10] = _crossHessian1[10][5] = VECTOR3(1,0,0);

  _crossHessian1[8][10] = _crossHessian1[10][8] =
  _crossHessian1[5][7] = _crossHessian1[7][5] =
  _crossHessian1[4][11] = _crossHessian1[11][4] = VECTOR3(-1,0,0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// cross product gradient of the left face
// so many zeros -- this is just begging to be refeactored into something more efficient
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 THETA_FASTER::crossHessian0(const int i, const int j) const
{
  return _crossHessian0[i][j];
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// cross product gradient of the left face
// so many zeros -- this is just begging to be refeactored into something more efficient
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 THETA_FASTER::crossHessian1(const int i, const int j) const
{
  return _crossHessian1[i][j];
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX12 THETA_FASTER::hessian(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  TIMER functionTimer("THETA_FASTER::hessian");
  
  // do cos theta
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const REAL cosTheta = n0.dot(n1);
  
  // do sine theta
  const VECTOR3& v1 = flap[1];
  const VECTOR3& v2 = flap[2];
  const VECTOR3 e = (v1 - v2) / (v1 - v2).norm();
  const REAL sinTheta = (n0.cross(n1)).dot(e);

  const MATRIX3x12 nGrad0 = normalGradient0(flap);
  const MATRIX3x12 nGrad1 = normalGradient1(flap);
  const MATRIX3x12 eGradient = edgeGradient(flap);

  // get the gradient of the cross-product
  MATRIX3x12 crossGradient;
  for (int x = 0; x < 12; x++)
  {
    const VECTOR3 left = nGrad0.col(x).cross(n1);
    const VECTOR3 right = n0.cross(nGrad1.col(x));
    crossGradient.col(x) = left + right;
  }
  
  // get the outer products first
  MATRIX12 cosThetaH = nGrad0.transpose() * nGrad1 + nGrad1.transpose() * nGrad0;
  const MATRIX12 CE = crossGradient.transpose() * eGradient;
  MATRIX12 sinThetaH = 0.5 * (CE + CE.transpose());

  const VECTOR3 e0 = flap[1] - flap[0];
  const VECTOR3 e1 = flap[2] - flap[0];
  const VECTOR3 e2 = flap[1] - flap[3];
  const VECTOR3 e3 = flap[2] - flap[3];
  const VECTOR3 z0 = e1.cross(e0);
  const VECTOR3 z1 = e2.cross(e3);

  const REAL invZnorm0 = 1.0 / z0.norm();
  const REAL invZnormCubed0 = invZnorm0 * invZnorm0 * invZnorm0;
  const REAL invZnormFive0  = invZnormCubed0 * invZnorm0 * invZnorm0;
  const REAL invZnorm1 = 1.0 / z1.norm();
  const REAL invZnormCubed1 = invZnorm1 * invZnorm1 * invZnorm1;
  const REAL invZnormFive1  = invZnormCubed1 * invZnorm1 * invZnorm1;

  const MATRIX3x12 zGradient0 = crossGradient0(flap);
  const MATRIX3x12 zGradient1 = crossGradient1(flap);

  // matrix is symmetric
  for (int j = 0; j < 12; j++)
    for (int i = j; i < 12; i++)
    {
      VECTOR3 h0ij;
      {
        const VECTOR3 zGradj0 = zGradient0.col(j);
        const VECTOR3 zGradi0 = zGradient0.col(i);
        const REAL zzi = z0.dot(zGradi0);
        const REAL zzj = z0.dot(zGradj0);
        const REAL alpha0 = zzi * invZnormCubed0;
        const VECTOR3 zHessian0 = crossHessian0(i,j);
        
        const REAL alphaGrad0 = invZnormCubed0 * ((zGradi0.dot(zGradj0)) + z0.dot(zHessian0)) - 
                                invZnormFive0 * 3.0 * zzi * zzj;
        h0ij = -1.0 * invZnormCubed0 * zzj * zGradi0 +
               invZnorm0 * zHessian0 - alpha0 * zGradj0 - alphaGrad0 * z0;
      }
      VECTOR3 h1ij;
      {
        const VECTOR3 zGradj1 = zGradient1.col(j);
        const VECTOR3 zGradi1 = zGradient1.col(i);
        const REAL zzi = z1.dot(zGradi1);
        const REAL zzj = z1.dot(zGradj1);
        const REAL alpha1 = zzi * invZnormCubed1;
        const VECTOR3 zHessian1 = crossHessian1(i,j);

        const REAL alphaGrad1 = invZnormCubed1 * (zGradi1.dot(zGradj1) + z1.dot(zHessian1)) - 
                                invZnormFive1 * 3.0 * zzi * zzj;
        h1ij = -1.0 * invZnormCubed1 * zzj * zGradi1 +
               invZnorm1 * zHessian1 - alpha1 * zGradj1 - alphaGrad1 * z1;
      }
      
      const VECTOR3 g0i = nGrad0.col(i);
      const VECTOR3 g0j = nGrad0.col(j);
      const VECTOR3 g1i = nGrad1.col(i);
      const VECTOR3 g1j = nGrad1.col(j);

      const VECTOR3 crossHessian = g0i.cross(g1j) + g0j.cross(g1i) + h0ij.cross(n1) + n0.cross(h1ij);
      const REAL sinTerm = crossHessian.dot(e);
      const REAL cosTerm = (h0ij.dot(n1) + h1ij.dot(n0));
      sinThetaH(i,j) += sinTerm;
      cosThetaH(i,j) += cosTerm;

      // don't symmetrize the diagonal
      if (i == j) continue;

      sinThetaH(j,i) += sinTerm;
      cosThetaH(j,i) += cosTerm;
    }

  return sinThetaH * cosTheta - cosThetaH * sinTheta;
}

} // SHELL
} // HOBAK
