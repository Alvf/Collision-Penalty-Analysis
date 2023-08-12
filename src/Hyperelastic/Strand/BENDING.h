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
// July 15, 2022 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef STRAND_BENDING_H
#define STRAND_BENDING_H

#include "SETTINGS.h"
#include "MATRIX_UTIL.h"

namespace HOBAK {
namespace STRAND {

///////////////////////////////////////////////////////////////////////////////////////////////////
// Abstract interface for this one is still a little ambiguous.
//
// The Bergou 2008 and 2010 models use a kappa term, but an isotropic model would just use a
// theta (bending angle) term. Going to veirfy against existing codes first before settling
// on a single interface.
///////////////////////////////////////////////////////////////////////////////////////////////////
class BENDING
{
public:
  BENDING(const REAL& mu) { _mu = mu; };
  virtual ~BENDING() { };

  // Computes the strain energy density
  virtual REAL psi(const VECTOR2& restKappa, const VECTOR2& kappa) const = 0;

  /*
  // f-based gradient
  virtual VECTOR3 PK1(const VECTOR3& f) const = 0;

  // f-based Hessian
  virtual MATRIX3 hessian(const VECTOR3& f) const = 0;
 
  // f-based clamped Hessian
  virtual MATRIX3 clampedHessian(const VECTOR3& f) const {
    return clampEigenvalues(hessian(f));
  };
  */
  virtual VECTOR9 spatialGradient(const VECTOR2& restKappa,
                                  const VECTOR2& kappa,
                                  const VECTOR3& tangent,
                                  const VECTOR3& previousTangent,
                                  const REAL& restLength,
                                  const REAL& previousRestLength,
                                  const REAL& length,
                                  const REAL& previousLength,
                                  const MATRIX2& perCornerKappa,
                                  const std::vector<VECTOR3>& director1,
                                  const std::vector<VECTOR3>& director2,
                                  const VECTOR3& kb) const = 0;

  // The name of the material
  virtual std::string name() const = 0;

protected:
  REAL _mu;
};

}
}

#endif
