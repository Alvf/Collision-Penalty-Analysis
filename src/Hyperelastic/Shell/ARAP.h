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

#ifndef SHELL_ARAP_H
#define SHELL_ARAP_H

#include "STRETCHING.h"

namespace HOBAK {
namespace SHELL {

class ARAP : public STRETCHING
{
public:
  ARAP(const REAL& mu, const REAL& lambda);

  // Computes the strain energy density
  virtual REAL psi(const MATRIX3x2& F) const override;

  // Computes the first Piola-Kirchoff PK1 stress
  virtual MATRIX3x2 PK1(const MATRIX3x2& F) const override;

  // Computes the derivative of the PK1 stress
  virtual MATRIX6 hessian(const MATRIX3x2& F) const override;

  // Computes the derivative of the PK1 stress, clamped to semi-positive definiteness
  virtual MATRIX6 clampedHessian(const MATRIX3x2& F) const override;
  
  // The name of the material
  virtual std::string name() const override;
};

}
}

#endif
