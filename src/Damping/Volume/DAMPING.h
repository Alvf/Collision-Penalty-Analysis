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
#ifndef VOLUME_DAMPING_H
#define VOLUME_DAMPING_H

#include "SETTINGS.h"

namespace HOBAK {
namespace VOLUME {

class DAMPING
{
public:
  // Computes the strain energy density
  virtual REAL psi(const MATRIX3& F, const MATRIX3& Fdot) const = 0;

  // Computes the first Piola-Kirchoff PK1 stress
  virtual MATRIX3 PK1(const MATRIX3& F, const MATRIX3& Fdot) const = 0;

  // The name of the material
  virtual std::string name() const = 0;

  // Computes the derivative of the PK1 stress
  virtual MATRIX9 hessian(const MATRIX3& F, const MATRIX3& Fdot) const = 0;

  // Computes the derivative of the PK1 stress, clamped to semi-positive definiteness
  virtual MATRIX9 clampedHessian(const MATRIX3& F, const MATRIX3& Fdot) const = 0;

  // The asymmetric term in damping that can occur because we take
  // a mixed derivative w.r.t. position and velocity. Still not sure how
  // valuable it is, but will include it here so we can do some testing
  virtual MATRIX9 positionGradient(const MATRIX3& F, const MATRIX3& Fdot) const { return MATRIX9::Zero(); };

  REAL& mu() { return _mu; };
  const REAL mu() const { return _mu; };

protected:
  REAL _mu;
};

}
}

#endif
