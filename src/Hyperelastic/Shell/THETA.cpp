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

#include "THETA.h"
#include "util/MATRIX_UTIL.h"
#include "TIMER.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace SHELL {

THETA::THETA(const REAL& mu) :
  BENDING(mu)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string THETA::name() const
{ 
  return std::string("THETA"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR12 flatten(const vector<VECTOR3>& flap)
{
  assert(flap.size() == 4);
  VECTOR12 flattened;
  int index = 0;
  for (int x = 0; x < 4; x++)
    for (int y = 0; y < 3; y++, index++)
      flattened[index] = flap[x][y];

  return flattened;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// normal of the left face
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 THETA::normal0(const vector<VECTOR3>& flap) const
{
  const VECTOR3 e0 = flap[1] - flap[0];
  const VECTOR3 e1 = flap[2] - flap[0];
  const VECTOR3 crossed = e1.cross(e0);
  return crossed / crossed.norm();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x12 THETA::edgeGradient(const vector<VECTOR3>& flap) const
{
  const VECTOR3& v1 = flap[1];
  const VECTOR3& v2 = flap[2];
  const VECTOR3 z = (v1 - v2);

  MATRIX3x12 zGradient;
  zGradient.setZero();

  zGradient.col(3) = VECTOR3(1,0,0); 
  zGradient.col(4) = VECTOR3(0,1,0); 
  zGradient.col(5) = VECTOR3(0,0,1);

  zGradient.col(6) = VECTOR3(-1,0,0); 
  zGradient.col(7) = VECTOR3(0,-1,0); 
  zGradient.col(8) = VECTOR3(0,0,-1);

  MATRIX3x12 gradient;
  const REAL invZnorm = 1.0 / z.norm();
  const REAL invZnormCubed = 1.0 / pow(z.dot(z), 1.5);
  for (int i = 0; i < 12; i++)
  {
    const VECTOR3 zg = zGradient.col(i);
    gradient.col(i) = invZnorm * zg - invZnormCubed * (z.dot(zg)) * z;
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// cross gradient of the right face
// so many zeros -- this is just begging to be refeactored into something more efficient
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x12 THETA::crossGradient1(const vector<VECTOR3>& flap)
{
  const VECTOR12 x = flatten(flap);
  MATRIX3x12 zGradient;
  zGradient.setZero();
  zGradient.col(3) = VECTOR3( 0,
                              x[11] - x[8],
                             -x[10] + x[7]);

  zGradient.col(4) = VECTOR3(-x[11] + x[8],
                              0,
                             -x[6] + x[9]);

  zGradient.col(5) = VECTOR3( x[10] - x[7],
                              x[6] - x[9],
                              0);

  zGradient.col(6) = VECTOR3( 0,
                             -x[11] + x[5],
                              x[10] - x[4]);

  zGradient.col(7) = VECTOR3( x[11] - x[5],
                              0,
                              x[3] - x[9]);

  zGradient.col(8) = VECTOR3(-x[10] + x[4],
                             -x[3] + x[9],
                              0);
  
  zGradient.col(9) = VECTOR3( 0,
                             -x[5] + x[8],
                              x[4] - x[7]);
  
  zGradient.col(10) = VECTOR3( x[5] - x[8],
                               0,
                              -x[3] + x[6]);

  zGradient.col(11) = VECTOR3(-x[4] + x[7],
                               x[3] - x[6],
                               0);

  return zGradient;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// normal gradient of the right face
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x12 THETA::normalGradient1(const vector<VECTOR3>& flap) const
{
  const VECTOR3 e2 = flap[1] - flap[3];
  const VECTOR3 e3 = flap[2] - flap[3];
  const VECTOR3 z = e2.cross(e3);

  const MATRIX3x12 zGradient = crossGradient1(flap);

  MATRIX3x12 gradient;
  const REAL invZnorm = 1.0 / z.norm();
  const REAL invZnormCubed = 1.0 / pow(z.dot(z), 1.5);
  for (int i = 0; i < 12; i++)
  {
    const VECTOR3 zg = zGradient.col(i);
    gradient.col(i) = invZnorm * zg - invZnormCubed * (z.dot(zg)) * z;
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// cross product gradient of the left face
// so many zeros -- this is just begging to be refeactored into something more efficient
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 THETA::crossHessian0(const int i, const int j) const
{
  if (((i == 0 && j == 4)) || ((i == 4) && (j == 0)) ||
      ((i == 1 && j == 6)) || ((i == 6) && (j == 1)) ||
      ((i == 3 && j == 7)) || ((i == 7) && (j == 3)))
    return VECTOR3(0,0,-1);
  
  if (((i == 0 && j == 5)) || ((i == 5) && (j == 0)) ||
      ((i == 2 && j == 6)) || ((i == 6) && (j == 2)) ||
      ((i == 3 && j == 8)) || ((i == 8) && (j == 3)))
    return VECTOR3(0,1,0);
  
  if (((i == 0 && j == 7)) || ((i == 7) && (j == 0)) ||
      ((i == 1 && j == 3)) || ((i == 3) && (j == 1)) ||
      ((i == 4 && j == 6)) || ((i == 6) && (j == 4)))
    return VECTOR3(0,0,1);
 
  if (((i == 0 && j == 8)) || ((i == 8) && (j == 0)) ||
      ((i == 2 && j == 3)) || ((i == 3) && (j == 2)) ||
      ((i == 5 && j == 6)) || ((i == 6) && (j == 5)))
    return VECTOR3(0,-1,0);

  if (((i == 1 && j == 5)) || ((i == 5) && (j == 1)) ||
      ((i == 2 && j == 7)) || ((i == 7) && (j == 2)) ||
      ((i == 4 && j == 8)) || ((i == 8) && (j == 4)))
    return VECTOR3(-1,0,0);
  
  if (((i == 1 && j == 8)) || ((i == 8) && (j == 1)) ||
      ((i == 2 && j == 4)) || ((i == 4) && (j == 2)) ||
      ((i == 5 && j == 7)) || ((i == 7) && (j == 5)))
    return VECTOR3(1,0,0);

  return VECTOR3(0,0,0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// cross product gradient of the left face
// so many zeros -- this is just begging to be refeactored into something more efficient
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 THETA::crossHessian1(const int i, const int j) const
{
  if (((i == 3 && j == 7)) || ((i == 7) && (j == 3)) ||
      ((i == 4 && j == 9)) || ((i == 9) && (j == 4)) ||
      ((i == 6 && j == 10)) || ((i == 10) && (j == 6)))
    return VECTOR3(0,0,1);

  if (((i == 3 && j == 8)) || ((i == 8) && (j == 3)) ||
      ((i == 6 && j == 11)) || ((i == 11) && (j == 6)) ||
      ((i == 5 && j == 9)) || ((i == 9) && (j == 5)))
    return VECTOR3(0,-1,0);

  if (((i == 7 && j == 9)) || ((i == 9) && (j == 7)) ||
      ((i == 4 && j == 6)) || ((i == 6) && (j == 4)) ||
      ((i == 3 && j == 10)) || ((i == 10) && (j == 3)))
    return VECTOR3(0,0,-1);

  if (((i == 8 && j == 9)) || ((i == 9) && (j == 8)) ||
      ((i == 5 && j == 6)) || ((i == 6) && (j == 5)) ||
      ((i == 3 && j == 11)) || ((i == 11) && (j == 3)))
    return VECTOR3(0,1,0);

  if (((i == 4 && j == 8)) || ((i == 8) && (j == 4)) ||
      ((i == 7 && j == 11)) || ((i == 11) && (j == 7)) ||
      ((i == 5 && j == 10)) || ((i == 10) && (j == 5)))
    return VECTOR3(1,0,0);

  if (((i == 8 && j == 10)) || ((i == 10) && (j == 8)) ||
      ((i == 5 && j == 7)) || ((i == 7) && (j == 5)) ||
      ((i == 4 && j == 11)) || ((i == 11) && (j == 4)))
    return VECTOR3(-1,0,0);
  
  return VECTOR3(0,0,0);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// cross gradient of the left face
// so many zeros -- this is just begging to be refeactored into something more efficient
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x12 THETA::crossGradient0(const vector<VECTOR3>& flap)
{
  const VECTOR12 x = flatten(flap);
  MATRIX3x12 zGradient;
  zGradient.setZero();
  zGradient.col(0) = VECTOR3( 0,
                              x[5] - x[8],
                             -x[4] + x[7]);          
  zGradient.col(1) = VECTOR3(-x[5] + x[8],
                              0,
                              x[3] - x[6]);
  zGradient.col(2) = VECTOR3( x[4] - x[7],
                             -x[3] + x[6],
                              0);
  zGradient.col(3) = VECTOR3( 0,
                             -x[2] + x[8],
                              x[1] - x[7]);
  zGradient.col(4) = VECTOR3( x[2] - x[8],
                              0,
                             -x[0] + x[6]);
  zGradient.col(5) = VECTOR3(-x[1] + x[7],
                              x[0] - x[6],
                              0);
  zGradient.col(6) = VECTOR3( 0,
                              x[2] - x[5],
                             -x[1] + x[4]);
  zGradient.col(7) = VECTOR3(-x[2] + x[5],
                              0,
                              x[0] - x[3]);
  zGradient.col(8) = VECTOR3( x[1] - x[4],
                             -x[0] + x[3],
                              0);

  return zGradient;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// normal gradient of the left face
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x12 THETA::normalGradient0(const vector<VECTOR3>& flap) const
{
  const VECTOR3 e0 = flap[1] - flap[0];
  const VECTOR3 e1 = flap[2] - flap[0];
  const VECTOR3 z = e1.cross(e0);

  const MATRIX3x12 zGradient = crossGradient0(flap);

  MATRIX3x12 gradient;
  const REAL invZnorm = 1.0 / z.norm();
  const REAL invZnormCubed = 1.0 / pow(z.dot(z), 1.5);
  for (int i = 0; i < 12; i++)
  {
    const VECTOR3 zg = zGradient.col(i);
    gradient.col(i) = invZnorm * zg - invZnormCubed * (z.dot(zg)) * z;
  }

  return gradient;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// normal of the right face
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR3 THETA::normal1(const vector<VECTOR3>& flap) const
{
  const VECTOR3 e2 = flap[1] - flap[3];
  const VECTOR3 e3 = flap[2] - flap[3];
  const VECTOR3 crossed = e2.cross(e3);
  return crossed / crossed.norm();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL THETA::psi(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  // debugging the sine and cosine terms
  //return cosThetaPsi(flap);
  //return sinThetaPsi(flap);

  // do cos theta
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const REAL cosTheta = n0.dot(n1);
  
  // do sine theta
  const VECTOR3& v1 = flap[1];
  const VECTOR3& v2 = flap[2];
  const VECTOR3 e = (v1 - v2) / (v1 - v2).norm();
  const REAL sinTheta = (n0.cross(n1)).dot(e);

  // get the robust theta
  return atan2(sinTheta, cosTheta);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL THETA::cosThetaPsi(const vector<VECTOR3>& flap) const
{
  // do cos theta first
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  return n0.dot(n1);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL THETA::sinThetaPsi(const vector<VECTOR3>& flap) const
{
  // do sin theta first
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const VECTOR3& v1 = flap[1];
  const VECTOR3& v2 = flap[2];
  const VECTOR3 e = (v1 - v2) / (v1 - v2).norm();
  return (n0.cross(n1)).dot(e);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR12 THETA::cosThetaGradient(const vector<VECTOR3>& flap) const
{
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const MATRIX3x12 nGrad0 = normalGradient0(flap);
  const MATRIX3x12 nGrad1 = normalGradient1(flap);
  return nGrad0.transpose() * n1 + nGrad1.transpose() * n0;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR12 THETA::sinThetaGradient(const vector<VECTOR3>& flap) const
{
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const VECTOR3& v1 = flap[1];
  const VECTOR3& v2 = flap[2];
  const VECTOR3 e = (v1 - v2) / (v1 - v2).norm();

  const MATRIX3x12 nGrad0 = normalGradient0(flap);
  const MATRIX3x12 nGrad1 = normalGradient1(flap);

  // get the gradient of the cross-product
  MATRIX3x12 crossGradient;
  for (int x = 0; x < 12; x++)
  {
    const VECTOR3 nGrad0i = nGrad0.col(x);
    const VECTOR3 nGrad1i = nGrad1.col(x);
    const VECTOR3 left  = nGrad0i.cross(n1);
    const VECTOR3 right = n0.cross(nGrad1i);
    crossGradient.col(x) = left + right;
  }

  const MATRIX3x12 eGradient = edgeGradient(flap);
  const VECTOR3 crossed = n0.cross(n1);
  return crossGradient.transpose() * e + eGradient.transpose() * crossed;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR12 THETA::gradient(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  //return cosThetaGradient(flap);
  //return sinThetaGradient(flap);

  // do cos theta
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const REAL cosTheta = n0.dot(n1);
  
  // do sine theta
  const VECTOR3& v1 = flap[1];
  const VECTOR3& v2 = flap[2];
  const VECTOR3 e = (v1 - v2) / (v1 - v2).norm();
  const REAL sinTheta = (n0.cross(n1)).dot(e);

  return sinThetaGradient(flap) * cosTheta - cosThetaGradient(flap) * sinTheta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
vector<MATRIX12> THETA::normalHessian0(const vector<VECTOR3>& flap) const
{
  //const VECTOR3 n0 = normal0(flap);
  //const VECTOR3 n1 = normal1(flap);
  const VECTOR3 e0 = flap[1] - flap[0];
  const VECTOR3 e1 = flap[2] - flap[0];
  const VECTOR3 z = e1.cross(e0);

  vector<MATRIX12> result(3);
  for (int i = 0; i < 3; i++)
    result[i].setZero();
  
  const REAL invZnorm = 1.0 / z.norm();
  const REAL invZnormCubed = 1.0 / pow(z.dot(z), 1.5);
  const REAL invZnormFive  = 1.0 / pow(z.dot(z), 2.5);

  const MATRIX3x12 zGradient = crossGradient0(flap);

  for (int j = 0; j < 12; j++)
    for (int i = 0; i < 12; i++)
    {
      const VECTOR3 zGradj = zGradient.col(j);
      const VECTOR3 zGradi = zGradient.col(i);
      const REAL alpha = (z.dot(zGradi)) * invZnormCubed;
      const VECTOR3 zHessian = crossHessian0(i,j);

      const REAL alphaGrad = invZnormCubed * (zGradi.dot(zGradj)) + 
                             invZnormCubed * (z.dot(zHessian)) - 
                             invZnormFive * 3.0 * (z.dot(zGradi)) * (zGradj.dot(z));
      const VECTOR3 entry = -1.0 * invZnormCubed * (zGradj.dot(z)) * zGradi +
                            invZnorm * zHessian - alpha * zGradj - alphaGrad * z;

      result[0](i,j) = entry[0];
      result[1](i,j) = entry[1];
      result[2](i,j) = entry[2];
    }

  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
vector<MATRIX12> THETA::normalHessian1(const vector<VECTOR3>& flap) const
{
  //const VECTOR3 n0 = normal0(flap);
  //const VECTOR3 n1 = normal1(flap);
  const VECTOR3 e2 = flap[1] - flap[3];
  const VECTOR3 e3 = flap[2] - flap[3];
  const VECTOR3 z = e2.cross(e3);

  vector<MATRIX12> result(3);
  for (int i = 0; i < 3; i++)
    result[i].setZero();
  
  const REAL invZnorm = 1.0 / z.norm();
  const REAL invZnormCubed = 1.0 / pow(z.dot(z), 1.5);
  const REAL invZnormFive  = 1.0 / pow(z.dot(z), 2.5);

  const MATRIX3x12 zGradient = crossGradient1(flap);

  for (int j = 0; j < 12; j++)
    for (int i = 0; i < 12; i++)
    {
      const VECTOR3 zGradj = zGradient.col(j);
      const VECTOR3 zGradi = zGradient.col(i);
      const REAL alpha = (z.dot(zGradi)) * invZnormCubed;
      const VECTOR3 zHessian = crossHessian1(i,j);

      const REAL alphaGrad = invZnormCubed * (zGradi.dot(zGradj)) + invZnormCubed * (z.dot(zHessian)) - invZnormFive * 3.0 * (z.dot(zGradi)) * (zGradj.dot(z));
      VECTOR3 entry = -1.0 * invZnormCubed * (zGradj.dot(z)) * zGradi +
                      invZnorm * zHessian - alpha * zGradj - alphaGrad * z;

      result[0](i,j) = entry[0];
      result[1](i,j) = entry[1];
      result[2](i,j) = entry[2];
    }

  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX12 THETA::cosThetaHessian(const vector<VECTOR3>& flap) const
{
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const MATRIX3x12 nGrad0 = normalGradient0(flap);
  const MATRIX3x12 nGrad1 = normalGradient1(flap);

  // get the outer products first
  MATRIX12 result = nGrad0.transpose() * nGrad1 + nGrad1.transpose() * nGrad0;

  const vector<MATRIX12> hessian0 = normalHessian0(flap);
  const vector<MATRIX12> hessian1 = normalHessian1(flap);

  for (int y = 0; y < 12; y++)
    for (int x = 0; x < 12; x++)
    {
      VECTOR3 H0(hessian0[0](x,y), hessian0[1](x,y), hessian0[2](x,y));
      VECTOR3 H1(hessian1[0](x,y), hessian1[1](x,y), hessian1[2](x,y));
      //result(x,y) += normalHessian0(flap, x,y).dot(n1) + normalHessian1(flap, x,y).dot(n0);
      result(x,y) += H0.dot(n1) + H1.dot(n0);
    }

  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX12 THETA::sinThetaHessian(const vector<VECTOR3>& flap) const
{
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const MATRIX3x12 nGrad0 = normalGradient0(flap);
  const MATRIX3x12 nGrad1 = normalGradient1(flap);

  // get the gradient of the cross-product
  MATRIX3x12 crossGradient;
  for (int x = 0; x < 12; x++)
  {
    const VECTOR3 left = nGrad0.col(x).cross(n1);
    const VECTOR3 right = n0.cross(nGrad1.col(x));
    crossGradient.col(x) = left + right;
  }

  const MATRIX3x12 eGradient = edgeGradient(flap);
  const vector<MATRIX12> hessian0 = normalHessian0(flap);
  const vector<MATRIX12> hessian1 = normalHessian1(flap);

  const VECTOR3& v1 = flap[1];
  const VECTOR3& v2 = flap[2];
  const VECTOR3 e = (v1 - v2) / (v1 - v2).norm();

  MATRIX12 result;
  result.setZero();

  for (int j = 0; j < 12; j++)
    for (int i = 0; i < 12; i++)
    {
      const VECTOR3 h0ij(hessian0[0](i,j), hessian0[1](i,j), hessian0[2](i,j));
      const VECTOR3 h1ij(hessian1[0](i,j), hessian1[1](i,j), hessian1[2](i,j));
      
      const VECTOR3 g0i = nGrad0.col(i);
      const VECTOR3 g0j = nGrad0.col(j);
      const VECTOR3 g1i = nGrad1.col(i);
      const VECTOR3 g1j = nGrad1.col(j);

      const VECTOR3 crossHessian = g0i.cross(g1j) + g0j.cross(g1i) + h0ij.cross(n1) + n0.cross(h1ij);
      result(i,j) += crossHessian.dot(e);
    }

  result += 0.5 * (crossGradient.transpose() * eGradient + eGradient.transpose() * crossGradient);
  return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX12 THETA::hessian(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  TIMER functionTimer("THETA::hessian");
  //return cosThetaHessian(flap);
  //return sinThetaHessian(flap);
  
  // do cos theta
  const VECTOR3 n0 = normal0(flap);
  const VECTOR3 n1 = normal1(flap);
  const REAL cosTheta = n0.dot(n1);
  
  // do sine theta
  const VECTOR3& v1 = flap[1];
  const VECTOR3& v2 = flap[2];
  const VECTOR3 e = (v1 - v2) / (v1 - v2).norm();
  const REAL sinTheta = (n0.cross(n1)).dot(e);

  return sinThetaHessian(flap) * cosTheta - cosThetaHessian(flap) * sinTheta;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Compute the bending angle based on a rest flap.
//
// In the case of a bending spring, it's just the difference between opposing vertices
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL THETA::restAngle(const std::vector<VECTOR3>& restFlap) const
{
  const VECTOR3& v0 = restFlap[0];
  const VECTOR3& v1 = restFlap[1];
  const VECTOR3& v2 = restFlap[2];
  const VECTOR3& v3 = restFlap[3];

  const VECTOR3 e20 = v2 - v0;
  const VECTOR3 e10 = v1 - v0;
  const VECTOR3 n0 = e20.cross(e10) / (e20 - e10).norm();

  const VECTOR3 e13 = v1 - v3;
  const VECTOR3 e23 = v2 - v3;
  const VECTOR3 n1 = e13.cross(e23) / (e13 - e23).norm();

  const VECTOR3 e12 = (v1 - v2) / (v1 - v2).norm();

  const REAL sinTheta = (n0.cross(n1)).dot(e12);
  const REAL cosTheta = n0.dot(n1);

  return atan2(sinTheta, cosTheta);
}

} // SHELL
} // HOBAK
