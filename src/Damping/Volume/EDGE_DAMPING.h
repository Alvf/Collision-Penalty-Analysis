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
#ifndef VOLUME_EDGE_DAMPING_H
#define VOLUME_EDGE_DAMPING_H

#include "SETTINGS.h"
#include <vector>

namespace HOBAK {
namespace VOLUME {

///////////////////////////////////////////////////////////////////////////////////
// This damping term doesn't actually seem to be written down anywhere,
// but is the obvious analog to the edge force described in
//
// "Dynamic deformables: implementation and production practicalities"
// by Kim and Eberle 2020. It complements the EDGE_COLLISION class, so if
// you want more details, go look at that header.
///////////////////////////////////////////////////////////////////////////////////
class EDGE_DAMPING
{
public:
  EDGE_DAMPING(const REAL& mu, const REAL& integratorConstant);
  ~EDGE_DAMPING() {};

  // get the strain energy
  virtual REAL psi(const std::vector<VECTOR3>& v,
                   const std::vector<VECTOR3>& vDot,
                   const VECTOR2& a, const VECTOR2& b) const;

  // This is the *gradient* of psi. The force is the *negative* gradient of psi.
  virtual VECTOR12 gradient(const std::vector<VECTOR3>& v,
                            const std::vector<VECTOR3>& vDot,
                            const VECTOR2& a, const VECTOR2& b) const;

  virtual std::string name() const;

  virtual MATRIX12 hessian(const std::vector<VECTOR3>& v,
                           const std::vector<VECTOR3>& vDot,
                           const VECTOR2& a,
                           const VECTOR2& b) const;

  virtual MATRIX12 clampedHessian(const std::vector<VECTOR3>& v,
                                  const std::vector<VECTOR3>& vDot,
                                  const VECTOR2& a,
                                  const VECTOR2& b) const;
    
protected:
  // convert the 12-vector in a way that imposes a consistent tet 
  // ordering for vertices and edges
  static void getVerticesAndEdges(const VECTOR12& x, 
                                  std::vector<VECTOR3>& v, 
                                  std::vector<VECTOR3>& e);

  // gradient of spring length, n' * (vDot[0] - vf)
  VECTOR12 normalVelocityGradient(const std::vector<VECTOR3>& vDot,
                                  const std::vector<VECTOR3>& e,
                                  const VECTOR3& n,
                                  const VECTOR2& a,
                                  const VECTOR2& b) const;

  // hessian of spring length, n' * (v[0] - v[2])
  MATRIX12 normalVelocityHessian(const std::vector<VECTOR3>& v,
                                 const std::vector<VECTOR3>& e,
                                 const VECTOR3& n,
                                 const VECTOR2& a,
                                 const VECTOR2& b) const;

  // partial of (vDot[0] - vf)
  static MATRIX3x12 vDiffPartial(const VECTOR2& a, const VECTOR2& b);

  // collision stiffness
  REAL _mu;

  // depending on whether you're using Newmark or BDF-1, a new constant
  // appears in front of \frac{\partial \dot\vv_- \dot\vv_f}{\partial \xx},
  // do it needs to be stored here
  REAL _integratorConstant;
};

} // VOLUME
} // HOBAK

#endif
