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
#ifndef VOLUME_EBERLE_DAMPING_H
#define VOLUME_EBERLE_DAMPING_H

#include "SETTINGS.h"
#include <vector>

namespace HOBAK {
namespace VOLUME {

///////////////////////////////////////////////////////////////////////////////////
// This is the damping term in the energy described from our SIGGRAPH 2020 notes
//
// "Dynamic deformables: implementation and production practicalities"
//
// by Kim and Eberle. This damping energy is specifically based on Dave's
// damping force, which is given in Eqn. 11.1. I don't think Dave is the one
// who came up with this term (he didn't claim to -- it's your usual dashpot), 
// and I don't think I've seen an energy form ever written down somewhere.
//
// The energy-based form here appears to be novel. Note that the force this
// produces is actually subtly different from the one Dave wrote down, so
// the force itself is slightly novel.
//
// Therefore, I get naming rights, and I'm calling this the Eberle Energy.
// tkim 8/27/21
///////////////////////////////////////////////////////////////////////////////////
class EBERLE_DAMPING
{
public:
  EBERLE_DAMPING(const REAL& mu, const REAL& integratorConstant);
  ~EBERLE_DAMPING() {};

  // get the strain energy
  virtual REAL psi(const std::vector<VECTOR3>& v,
                   const std::vector<VECTOR3>& vDot,
                   const VECTOR3& bary) const;

  // This is the *gradient* of psi. The force is the *negative* gradient of psi.
  virtual VECTOR12 gradient(const std::vector<VECTOR3>& v,
                            const std::vector<VECTOR3>& vDot,
                            const VECTOR3& bary) const;

  virtual std::string name() const;

  virtual MATRIX12 hessian(const std::vector<VECTOR3>& v,
                           const std::vector<VECTOR3>& vDot,
                           const VECTOR3& bary) const;

  virtual MATRIX12 clampedHessian(const std::vector<VECTOR3>& v,
                                  const std::vector<VECTOR3>& vDot,
                                  const VECTOR3& bary) const;
   
  virtual MATRIX12 velocityHessian(const std::vector<VECTOR3>& v,
                                   const std::vector<VECTOR3>& vDot,
                                   const VECTOR3& bary) const;

protected:
  // convert the 12-vector in a way that imposes a consistent tet 
  // ordering for vertices and edges
  static void getVerticesAndEdges(const VECTOR12& x, 
                                  std::vector<VECTOR3>& v, 
                                  std::vector<VECTOR3>& e);

  // get edges from vertices using a consistent ordering
  static void getEdges(const std::vector<VECTOR3>& v, 
                       std::vector<VECTOR3>& e);

  // gradient of spring length, n' * (vDot[0] - vf)
  VECTOR12 normalVelocityGradient(const std::vector<VECTOR3>& vDot,
                                  const std::vector<VECTOR3>& e,
                                  const VECTOR3& n,
                                  const VECTOR3& bary,
                                  const VECTOR3& vf) const;

  // hessian of spring length, n' * (v[0] - v[2])
  MATRIX12 normalVelocityHessian(const std::vector<VECTOR3>& v,
                                 const std::vector<VECTOR3>& e,
                                 const VECTOR3& n,
                                 const VECTOR3& bary) const;

  // partial of (vDot[0] - vf)
  static MATRIX3x12 vDiffPartial(const VECTOR3& bary);

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
