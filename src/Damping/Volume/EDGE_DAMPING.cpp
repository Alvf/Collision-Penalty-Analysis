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
#include "EDGE_DAMPING.h"
#include "MATRIX_UTIL.h"
#include "Collision/COLLISION_UTIL.h"
#include <iostream>

using namespace std;

namespace HOBAK {
namespace VOLUME {

EDGE_DAMPING::EDGE_DAMPING(const REAL& mu, const REAL& integratorConstant) :
    _mu(mu), _integratorConstant(integratorConstant)
{
}

std::string EDGE_DAMPING::name() const
{
  return "EDGE_DAMPING";
}

///////////////////////////////////////////////////////////////////////
// convert the 12-vector in a way that imposes a consistent tet 
// ordering for vertices and edges
///////////////////////////////////////////////////////////////////////
void EDGE_DAMPING::getVerticesAndEdges(const VECTOR12& x,
                                       vector<VECTOR3>& v,
                                       vector<VECTOR3>& e)
{
  v.resize(4);
  for (int i = 0; i < 4; i++)
  {
    v[i][0] = x[i * 3];
    v[i][1] = x[i * 3 + 1];
    v[i][2] = x[i * 3 + 2];
  }

  e.resize(2);
  e[0] = v[1] - v[0];
  e[1] = v[3] - v[2];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL EDGE_DAMPING::psi(const vector<VECTOR3>& v,
                       const vector<VECTOR3>& vDot,
                       const VECTOR2& a,
                       const VECTOR2& b) const
{
  // convert to edges
  vector<VECTOR3> e(2);
  e[0] = v[1] - v[0];
  e[1] = v[3] - v[2];

  // get the normal
  VECTOR3 n = e[1].cross(e[0]);
  n = n / n.norm();

  // get the interpolated vertices
  const VECTOR3 vaDot = (a[0] * vDot[0] + a[1] * vDot[1]);
  const VECTOR3 vbDot = (b[0] * vDot[2] + b[1] * vDot[3]);
  const VECTOR3 diff = vaDot - vbDot;

  // get the magnitude of the velocity
  REAL normalVelocityMagnitude = (n.dot(diff));
  return _mu * normalVelocityMagnitude * normalVelocityMagnitude;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR12 EDGE_DAMPING::gradient(const vector<VECTOR3>& v,
                                const vector<VECTOR3>& vDot,
                                const VECTOR2& a,
                                const VECTOR2& b) const
{
  // convert to edges
  vector<VECTOR3> e(2);
  e[0] = v[1] - v[0];
  e[1] = v[3] - v[2];

  // get the normal
  VECTOR3 n = e[1].cross(e[0]);
  n = n / n.norm();
 
  // get the interpolated vertices
  const VECTOR3 vaDot = (a[0] * vDot[0] + a[1] * vDot[1]);
  const VECTOR3 vbDot = (b[0] * vDot[2] + b[1] * vDot[3]);
  const VECTOR3 diff = vaDot - vbDot;

  // get the spring length, non-zero rest-length
  const REAL normalVelocityMagnitude = n.dot(diff);
  
  return 2.0 * _mu * normalVelocityMagnitude * normalVelocityGradient(vDot,e,n,a,b);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 EDGE_DAMPING::clampedHessian(const vector<VECTOR3>& v,
                                      const vector<VECTOR3>& vDot,
                                      const VECTOR2& a,
                                      const VECTOR2& b) const
{
  const MATRIX12 H = hessian(v,vDot,a,b);
  //return clampEigenvaluesToSemiNegative(H);
  return clampEigenvalues(H);
}

///////////////////////////////////////////////////////////////////////
// partial of (va - vb)
///////////////////////////////////////////////////////////////////////
MATRIX3x12 EDGE_DAMPING::vDiffPartial(const VECTOR2& a, const VECTOR2& b)
{
  MATRIX3x12 tPartial;
  tPartial.setZero();
  tPartial(0,0) = tPartial(1,1)  = tPartial(2,2) = a[0];
  tPartial(0,3) = tPartial(1,4)  = tPartial(2,5) = a[1];
  tPartial(0,6) = tPartial(1,7)  = tPartial(2,8) = -b[0];
  tPartial(0,9) = tPartial(1,10) = tPartial(2,11) = -b[1];

  return tPartial;
}

///////////////////////////////////////////////////////////////////////
// gradient of normal velocity, n' * (vDot[0] - vDot_f)
///////////////////////////////////////////////////////////////////////
VECTOR12 EDGE_DAMPING::normalVelocityGradient(const std::vector<VECTOR3>& vDot,
                                              const std::vector<VECTOR3>& e,
                                              const VECTOR3& n,
                                              const VECTOR2& a,
                                              const VECTOR2& b) const
{
  //nPartial = normal_gradient(x);
  MATRIX3x12 nPartial = normalGradientEE(e);
  MATRIX3x12 vPartial = _integratorConstant * vDiffPartial(a,b);

  // get the interpolated vertices
  const VECTOR3 vaDot = (a[0] * vDot[0] + a[1] * vDot[1]);
  const VECTOR3 vbDot = (b[0] * vDot[2] + b[1] * vDot[3]);
  const VECTOR3 diff = vaDot - vbDot;

  return nPartial.transpose() * diff + vPartial.transpose() * n;
}

///////////////////////////////////////////////////////////////////////
// hessian of spring length, n' * (v[2] - v[0])
///////////////////////////////////////////////////////////////////////
MATRIX12 EDGE_DAMPING::normalVelocityHessian(const std::vector<VECTOR3>& vDot,
                                             const std::vector<VECTOR3>& e,
                                             const VECTOR3& n,
                                             const VECTOR2& a,
                                             const VECTOR2& b) const
{
  // get the interpolated vertices
  const VECTOR3 vaDot = (a[0] * vDot[0] + a[1] * vDot[1]);
  const VECTOR3 vbDot = (b[0] * vDot[2] + b[1] * vDot[3]);
  const VECTOR3 t = vaDot - vbDot;

  MATRIX3x12 tPartial = vDiffPartial(a,b);

  // remember to incorporate the Newmark or BDF constant
  tPartial *= _integratorConstant;

  //% mode-3 contraction
  //[nx ny nz] = normal_hessian(x);
  //final = nx * delta(1) + ny * delta(2) + nz * delta(3);
  vector<MATRIX12> normalH = normalHessianEE(e);

  MATRIX12 contracted = t[0] * normalH[0] + 
                        t[1] * normalH[1] + 
                        t[2] * normalH[2];
  
  //nGrad= normal_gradient(x);
  MATRIX3x12 nGrad = normalGradientEE(e);

  //product = nGrad' * vGrad;
  //final = final + product + product';
  MATRIX12 product = nGrad.transpose() * tPartial;

  return contracted + product + product.transpose();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 EDGE_DAMPING::hessian(const std::vector<VECTOR3>& v,
                               const std::vector<VECTOR3>& vDot,
                               const VECTOR2& a,
                               const VECTOR2& b) const
{
  // convert to edges
  vector<VECTOR3> e(2);
  e[0] = v[1] - v[0];
  e[1] = v[3] - v[2];
  
  // get the normal
  VECTOR3 n = e[1].cross(e[0]);
  n = n / n.norm();

  // get the interpolated vertices
  const VECTOR3 vaDot = (a[0] * vDot[0] + a[1] * vDot[1]);
  const VECTOR3 vbDot = (b[0] * vDot[2] + b[1] * vDot[3]);
  const VECTOR3 diff = vaDot - vbDot;
  //const VECTOR3 diff = vbDot - vaDot;

  // get the spring length, non-zero rest-length
  REAL normalVelocityMagnitude = n.dot(diff);

  // ndotGrad    = ndot_gradient(x);
  VECTOR12 normalVelocityGrad = normalVelocityGradient(vDot,e,n,a,b);

  // ndotHessian = ndot_hessian(x);
  MATRIX12 normalVelocityH = normalVelocityHessian(vDot,e,n,a,b);
 
  // final = 2 * k * (ndotGrad * ndotGrad' + ndot * ndotHessian);
  return 2.0 * _mu * (normalVelocityGrad * normalVelocityGrad.transpose() + 
                      normalVelocityMagnitude * normalVelocityH);
}

} // VOLUME
} // HOBAK
