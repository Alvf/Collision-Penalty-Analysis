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

#include "DIHEDRAL.h"
#include "util/MATRIX_UTIL.h"
#include "TIMER.h"

namespace HOBAK {
namespace SHELL {

using namespace std;

DIHEDRAL::DIHEDRAL(const REAL& mu) :
  THETA_FASTER(mu)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string DIHEDRAL::name() const
{ 
  return std::string("DIHEDRAL"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL DIHEDRAL::psi(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  const VECTOR3 p[4] = {flap[0], flap[1], flap[2], flap[3]};

  const VECTOR3 p01 = p[1] - p[0];
  const VECTOR3 p20 = p[0] - p[2];
  const VECTOR3 p12 = p[2] - p[1];
  const VECTOR3 p32 = p[2] - p[3];
  const VECTOR3 p13 = p[3] - p[1];

  const VECTOR3 n0hat = (p[1] - p[0]).cross(p[2] - p[0]);
  const VECTOR3 n1hat = (p[2] - p[3]).cross(p[1] - p[3]);
  const VECTOR3 ehat = p12.normalized();

  const REAL sinTheta = (n0hat.cross(n1hat)).dot(ehat);
  const REAL cosTheta = n0hat.dot(n1hat); 

  const REAL theta = atan2(sinTheta, cosTheta); 
  const REAL thetaRest = restTheta;

  return _mu * 0.5 * (theta - thetaRest) * (theta - thetaRest);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
VECTOR12 DIHEDRAL::gradient(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  const VECTOR3 p[4] = {flap[0], flap[1], flap[2], flap[3]};

  const VECTOR3 p01 = p[1] - p[0];
  const VECTOR3 p20 = p[0] - p[2];
  const VECTOR3 p12 = p[2] - p[1];
  const VECTOR3 p32 = p[2] - p[3];
  const VECTOR3 p13 = p[3] - p[1];

  const VECTOR3 n0hat = (p[1] - p[0]).cross(p[2] - p[0]);
  const VECTOR3 n1hat = (p[2] - p[3]).cross(p[1] - p[3]);
  const VECTOR3 ehat = p12.normalized();

  const REAL sinTheta = (n0hat.cross(n1hat)).dot(ehat);
  const REAL cosTheta = n0hat.dot(n1hat); 

  const REAL theta = atan2(sinTheta, cosTheta); 

  return -_mu * (theta - restTheta) * THETA::gradient(flap, restTheta);
}
// David Eberle-style, from the course notes
/*
{
  const VECTOR3 p[4] = {flap[0], flap[1], flap[2], flap[3]};

  const VECTOR3 p01 = p[1] - p[0];
  const VECTOR3 p20 = p[0] - p[2];
  const VECTOR3 p12 = p[2] - p[1];
  const VECTOR3 p32 = p[2] - p[3];
  const VECTOR3 p13 = p[3] - p[1];

  const VECTOR3 n0hat = (p[1] - p[0]).cross(p[2] - p[0]);
  const VECTOR3 n1hat = (p[2] - p[3]).cross(p[1] - p[3]);
  const VECTOR3 ehat = p12.normalized();

  const REAL sinTheta = (n0hat.cross(n1hat)).dot(ehat);
  const REAL cosTheta = n0hat.dot(n1hat); 

  const REAL invN0Mag = 1.0 / n0hat.norm();
  const REAL invN1Mag = 1.0 / n1hat.norm();

  const VECTOR3 n0hatXehat = n0hat.cross(ehat);
  const VECTOR3 ehatXn1hat = ehat.cross(n1hat);

  VECTOR3 dSinThetadP[4];
  dSinThetadP[0] = invN0Mag * p12.cross(ehatXn1hat);
  dSinThetadP[1] = invN1Mag * p32.cross(n0hatXehat) +
                   invN0Mag * p20.cross(ehatXn1hat);
  dSinThetadP[2] = invN1Mag * p13.cross(n0hatXehat) +
                   invN0Mag * p01.cross(ehatXn1hat);
  dSinThetadP[3] = invN1Mag * n0hatXehat.cross(p12);

  const VECTOR3 projectN0OutFromN1 = (MATRIX3::Identity() - n0hat * n0hat.transpose()) * n1hat;
  const VECTOR3 projectN1OutFromN0 = (MATRIX3::Identity() - n1hat * n1hat.transpose()) * n0hat;
 
  VECTOR3 dCosThetadP[4];
  dCosThetadP[0] = -invN0Mag * p12.cross(projectN0OutFromN1);
  dCosThetadP[1] = -invN0Mag * p20.cross(projectN0OutFromN1);
  dCosThetadP[1] += -invN1Mag * p32.cross(projectN1OutFromN0);
  dCosThetadP[2] = -invN0Mag * p01.cross(projectN0OutFromN1);
  dCosThetadP[2] += -invN1Mag * p13.cross(projectN1OutFromN0);
  dCosThetadP[3] = -invN1Mag * projectN1OutFromN0.cross(p12);

  VECTOR3 dThetadP[4];
  REAL thetaDot = 0.0;
  for (int i = 0; i < 4; i++)
  {
      dThetadP[i] = dSinThetadP[i];
      dThetadP[i] *= cosTheta;
      dThetadP[i] -= sinTheta * dCosThetadP[i];
  }

  const REAL theta = atan2(sinTheta, cosTheta); 
  const REAL thetaRest = restTheta;

  // Note: thetaDist does not account for windings 
  const REAL thetaDist = theta - thetaRest;
  //const REAL mag = -ksTilde * thetaDist + -kdTilde * thetaDot;
  //const REAL mag = -ksTilde * thetaDist; 
  const REAL mag = _mu * thetaDist; 
 
  VECTOR3 spatialK_dThetadPi;
  VECTOR3 kd_dThetadPi;
  int blockIdxII, blockIdxIJ;

  VECTOR3 force[4];
  for (int i = 0; i < 4; i++)
    force[i].setZero();

  MATRIX3 dfdx[4][4];
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      dfdx[i][j].setZero();

  for (int i = 0; i < 4; i++)
    force[i] += mag * dThetadP[i];

  VECTOR12 gradient;
  int index = 0;
  for (int j = 0; j < 4; j++)
    for (int i = 0; i < 3; i++, index++)
      gradient[index] = force[j][i];

  return gradient;
}
*/

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX12 DIHEDRAL::hessian(const vector<VECTOR3>& flap, const REAL& restTheta) const
{
  TIMER functionTimer("DIHEDRAL::hessian");
  const VECTOR3 p[4] = {flap[0], flap[1], flap[2], flap[3]};

  const VECTOR3 n0hat = (p[1] - p[0]).cross(p[2] - p[0]);
  const VECTOR3 n1hat = (p[2] - p[3]).cross(p[1] - p[3]);
  const VECTOR3 ehat = (p[2] - p[1]).normalized();

  const REAL sinTheta = (n0hat.cross(n1hat)).dot(ehat);
  const REAL cosTheta = n0hat.dot(n1hat); 

  const REAL theta = atan2(sinTheta, cosTheta); 

  const VECTOR12 g = THETA::gradient(flap, restTheta);
  const MATRIX12 H = -_mu * ((theta - restTheta) * THETA_FASTER::hessian(flap, restTheta) - g * g.transpose());
  return H;
  //return clampEigenvalues(H);
  //return -_mu * ((theta - restTheta) * THETA_FASTER::hessian(flap, restTheta) - g * g.transpose());
}
// David Eberle-style, Gauss-Newton approximation
/*
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
*/

///////////////////////////////////////////////////////////////////////////////////////////////////
// Compute the bending angle based on a rest flap.
//
// In the case of a bending spring, it's just the difference between opposing vertices
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL DIHEDRAL::restAngle(const std::vector<VECTOR3>& restFlap) const
{
  //return (restFlap[0] - restFlap[3]).norm();
  return THETA::restAngle(restFlap);
}

} // SHELL
} // HOBAK
