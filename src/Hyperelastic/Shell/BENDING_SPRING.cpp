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

#include "BENDING_SPRING.h"
#include "util/MATRIX_UTIL.h"

#include <iostream>
using namespace std;

namespace HOBAK {
namespace SHELL {

BENDING_SPRING::BENDING_SPRING(const REAL& mu) :
  BENDING(mu)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string BENDING_SPRING::name() const
{ 
  return std::string("BENDING_SPRING"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL BENDING_SPRING::psi(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  const VECTOR3& p0 = flap[0];
  const VECTOR3& p1 = flap[3];

  const VECTOR3 d = p1 - p0;
  const REAL dNorm = d.norm();

  return -_mu * 0.5 * (dNorm - restTheta) * (dNorm - restTheta);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR12 BENDING_SPRING::gradient(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  const VECTOR3& p0 = flap[0];
  const VECTOR3& p1 = flap[3];
  const VECTOR3 d = p1 - p0;
  const REAL dNorm = d.norm();

  VECTOR12 g;
  g.setZero();
  g.segment<3>(0) = -d;
  g.segment<3>(9) = d;

  return -(_mu * (dNorm - restTheta) / dNorm) * g;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX12 BENDING_SPRING::hessian(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  const VECTOR3& p0 = flap[0];
  const VECTOR3& p1 = flap[3];
  const VECTOR3 d = p1 - p0;
  const REAL dNorm = d.norm();

  VECTOR12 g;
  g.setZero();
  g.segment<3>(0) = -d;
  g.segment<3>(9) = d;
  g *= 1.0 / dNorm;

  MATRIX12 Z;
  Z.setZero();
  Z.block<3,3>(0,0) = MATRIX3::Identity();
  Z.block<3,3>(9,9) = MATRIX3::Identity();
  Z.block<3,3>(0,9) = -1.0 * MATRIX3::Identity();
  Z.block<3,3>(9,0) = -1.0 * MATRIX3::Identity();

  const MATRIX12 H = (1.0 / dNorm) * (Z - g * g.transpose());
  return -_mu * ((dNorm - restTheta) * H + (g * g.transpose()));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Compute the bending angle based on a rest flap.
//
// In the case of a bending spring, it's just the difference between opposing vertices
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL BENDING_SPRING::restAngle(const std::vector<VECTOR3>& restFlap) const
{
  return (restFlap[0] - restFlap[3]).norm();
}

} // SHELL
} // HOBAK
