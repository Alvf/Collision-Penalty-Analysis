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
#include "EBERLE_DAMPING.h"
#include "MATRIX_UTIL.h"
#include "Collision/COLLISION_UTIL.h"
#include <iostream>

#define NEGATE 1
#define VERBOSE 0
#define TRUNCATED_FORCE 0

using namespace std;

namespace HOBAK {
namespace VOLUME {

EBERLE_DAMPING::EBERLE_DAMPING(const REAL& mu, const REAL& integratorConstant) :
    _mu(mu), _integratorConstant(integratorConstant)
{
}

std::string EBERLE_DAMPING::name() const
{
  return "EBERLE_DAMPING";
}

///////////////////////////////////////////////////////////////////////
// get edges from vertices using a consistent ordering
///////////////////////////////////////////////////////////////////////
void EBERLE_DAMPING::getEdges(const std::vector<VECTOR3>& v, 
                              std::vector<VECTOR3>& e)
{
  e.resize(3);
  e[0] = v[3] - v[2];
  e[1] = v[0] - v[2];
  e[2] = v[1] - v[2];
}

///////////////////////////////////////////////////////////////////////
// convert the 12-vector in a way that imposes a consistent tet 
// ordering for vertices and edges
//
// OLD: assumes v2, x[6,7,8] is the collision vertex
// assumes v0, x[0,1,2] is the collision vertex
///////////////////////////////////////////////////////////////////////
void EBERLE_DAMPING::getVerticesAndEdges(const VECTOR12& x,
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

  e.resize(3);
  e[0] = v[3] - v[2];
  e[1] = v[0] - v[2];
  e[2] = v[1] - v[2];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL EBERLE_DAMPING::psi(const vector<VECTOR3>& v,
                         const vector<VECTOR3>& vDot,
                         const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> e;
  getEdges(v, e);
  //const VECTOR3 bary = getBarycentricCoordinates(v);

  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();

  // face velocity
  const VECTOR3 vf = bary[0] * vDot[1] + bary[1] * vDot[2] + bary[2] * vDot[3];

  // get the magnitude of the velocity
#if NEGATE
  REAL normalVelocityMagnitude = (n.dot(vf - vDot[0]));
#else
  REAL normalVelocityMagnitude = (n.dot(vDot[0] - vf));
#endif
  return _mu * normalVelocityMagnitude * normalVelocityMagnitude;
}

///////////////////////////////////////////////////////////////////////
// This is the *gradient* of psi. The force is the *negative* gradient
// of psi.
///////////////////////////////////////////////////////////////////////
VECTOR12 EBERLE_DAMPING::gradient(const vector<VECTOR3>& v,
                                  const vector<VECTOR3>& vDot,
                                  const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> e;
  getEdges(v, e);
  
  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();
  
  // face velocity
  //const VECTOR3 bary = getBarycentricCoordinates(v);
  const VECTOR3 vf = bary[0] * vDot[1] + bary[1] * vDot[2] + bary[2] * vDot[3];
  
  // get the spring length, non-zero rest-length
#if NEGATE
  REAL normalVelocityMagnitude = (n.dot(vf - vDot[0]));
#else
  REAL normalVelocityMagnitude = (n.dot(vDot[0] - vf));
#endif

#if VERBOSE
  VECTOR3 diff = vDot[0] - vf;
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " diff:      " << diff.transpose() << endl;
  cout << " constant:  " << _integratorConstant << endl;
  cout << " magnitude: " << normalVelocityMagnitude << endl;
  cout << " gradient:  " << normalVelocityGradient(vDot,e,n,bary,vf) << endl;
#endif

#if TRUNCATED_FORCE
  return 2.0 * _mu * normalVelocityMagnitude * vDiffPartial(bary).transpose() * n;
#else
  return 2.0 * _mu * normalVelocityMagnitude * normalVelocityGradient(vDot,e,n,bary,vf);
#endif
}

///////////////////////////////////////////////////////////////////////
// This is the *gradient* of psi. The force is the *negative* gradient
// of psi.
///////////////////////////////////////////////////////////////////////
MATRIX12 EBERLE_DAMPING::velocityHessian(const vector<VECTOR3>& v,
                                         const vector<VECTOR3>& vDot,
                                         const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> e;
  getEdges(v, e);
  
  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();
  
  // face velocity
  //const VECTOR3 bary = getBarycentricCoordinates(v);
  const VECTOR3 vf = bary[0] * vDot[1] + bary[1] * vDot[2] + bary[2] * vDot[3];
  
  const VECTOR12 nGradient = normalVelocityGradient(vDot,e,n,bary,vf);
  return 2.0 * _mu * (nGradient * nGradient.transpose());
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 EBERLE_DAMPING::clampedHessian(const vector<VECTOR3>& v,
                                        const vector<VECTOR3>& vDot,
                                        const VECTOR3& bary) const
{
  return hessian(v, vDot, bary);
  //return clampEigenvalues(hessian(v, vDot, bary));
  //return clampEigenvalues(hessian(v, vDot));
  //return clampEigenvaluesToSemiNegative(hessian(v, vDot));
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX3x12 EBERLE_DAMPING::vDiffPartial(const VECTOR3& bary)
{
  MATRIX3x12 tPartial;
  tPartial.setZero();
#if NEGATE
  tPartial(0,0) = tPartial(1,1)  = tPartial(2,2) = -1.0;
  tPartial(0,3) = tPartial(1,4)  = tPartial(2,5) = bary[0];
  tPartial(0,6) = tPartial(1,7)  = tPartial(2,8) = bary[1];
  tPartial(0,9) = tPartial(1,10) = tPartial(2,11) = bary[2];
#else
  tPartial(0,0) = tPartial(1,1)  = tPartial(2,2) = 1.0;
  tPartial(0,3) = tPartial(1,4)  = tPartial(2,5) = -bary[0];
  tPartial(0,6) = tPartial(1,7)  = tPartial(2,8) = -bary[1];
  tPartial(0,9) = tPartial(1,10) = tPartial(2,11) = -bary[2];
#endif

  return tPartial;
}

///////////////////////////////////////////////////////////////////////
// gradient of normal velocity, n' * (vDot[0] - vDot_f)
///////////////////////////////////////////////////////////////////////
VECTOR12 EBERLE_DAMPING::normalVelocityGradient(const std::vector<VECTOR3>& vDot,
                                                const std::vector<VECTOR3>& e,
                                                const VECTOR3& n,
                                                const VECTOR3& bary,
                                                const VECTOR3& vf) const
{
  //nPartial = normal_gradient(x);
  MATRIX3x12 nPartial = normalGradientVF(e);
  MATRIX3x12 vPartial = _integratorConstant * vDiffPartial(bary);

  //f = nPartial' * (v2 - v0) + vPartial' * n;
  //return nPartial.transpose() * (v[2] - v[0]) + vPartial.transpose() * n;
#if NEGATE
#if VERBOSE
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  cout << " term1: " << (nPartial.transpose() * (vDot[0] - vf)).transpose() << endl;
  cout << " term2: " << (vPartial.transpose() * n).transpose() << endl;
#endif
  return nPartial.transpose() * (vf - vDot[0]) + vPartial.transpose() * n;
#else
  return nPartial.transpose() * (vDot[0] - vf) + vPartial.transpose() * n;
#endif
}

///////////////////////////////////////////////////////////////////////
// hessian of spring length, n' * (v[2] - v[0])
///////////////////////////////////////////////////////////////////////
MATRIX12 EBERLE_DAMPING::normalVelocityHessian(const std::vector<VECTOR3>& vDot,
                                               const std::vector<VECTOR3>& e,
                                               const VECTOR3& n,
                                               const VECTOR3& bary) const
{
  // remember we had to reorder vertices in a wonky way
  const VECTOR3 vf = bary[0] * vDot[1] + bary[1] * vDot[2] + bary[2] * vDot[3];
#if NEGATE
  const VECTOR3 t = vf - vDot[0];
#else
  const VECTOR3 t = vDot[0] - vf;
#endif

  MATRIX3x12 tPartial = vDiffPartial(bary);

  // remember to incorporate the Newmark or BDF constant
  tPartial *= _integratorConstant;

  //% mode-3 contraction
  //[nx ny nz] = normal_hessian(x);
  //final = nx * delta(1) + ny * delta(2) + nz * delta(3);
  vector<MATRIX12> normalH = normalHessianVF(e);

  MATRIX12 contracted = t[0] * normalH[0] + 
                        t[1] * normalH[1] + 
                        t[2] * normalH[2];
  
  //nGrad= normal_gradient(x);
  MATRIX3x12 nGrad = normalGradientVF(e);

  //product = nGrad' * vGrad;
  //final = final + product + product';
  MATRIX12 product = nGrad.transpose() * tPartial;

  return contracted + product + product.transpose();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX12 EBERLE_DAMPING::hessian(const std::vector<VECTOR3>& v,
                                 const std::vector<VECTOR3>& vDot,
                                 const VECTOR3& bary) const
{
  // convert to vertices and edges
  vector<VECTOR3> e;
  getEdges(v,e);
  
  // get the normal
  VECTOR3 n = e[2].cross(e[0]);
  n = n / n.norm();

  // face velocity
  //const VECTOR3 bary = getBarycentricCoordinates(v);
  const VECTOR3 vf = bary[0] * vDot[1] + bary[1] * vDot[2] + bary[2] * vDot[3];

  // get the spring length, non-zero rest-length
#if NEGATE
  REAL normalVelocityMagnitude = n.dot(vf - vDot[0]);
#else
  REAL normalVelocityMagnitude = n.dot(vDot[0] - vf);
#endif

  // ndotGrad    = ndot_gradient(x);
  VECTOR12 normalVelocityGrad = normalVelocityGradient(vDot,e,n,bary,vf);

  // ndotHessian = ndot_hessian(x);
  MATRIX12 normalVelocityH = normalVelocityHessian(vDot,e,n,bary);

#if TRUNCATED_FORCE
  MATRIX3x12 vPartial = vDiffPartial(bary);
  MATRIX3x12 nPartial = normalGradientVF(e);
  // final = 2 * k * (ndotGrad * ndotGrad' + ndot * ndotHessian);
  //
  VECTOR12 force = vPartial.transpose() * n;

  MATRIX12 product = vPartial.transpose() * nPartial;

  //cout << " first: " << endl << normalVelocityGrad * force.transpose() << endl;
  //cout << " second : " << endl << normalVelocityMagnitude * product << endl;

  return 2.0 * _mu * (normalVelocityGrad * force.transpose() +
                      normalVelocityMagnitude * product);
                      //normalVelocityMagnitude * (product + product.transpose()));
#else
  // final = 2 * k * (ndotGrad * ndotGrad' + ndot * ndotHessian);
  return 2.0 * _mu * (normalVelocityGrad * normalVelocityGrad.transpose() + 
                      normalVelocityMagnitude * normalVelocityH);
#endif
}

} // VOLUME
} // HOBAK
