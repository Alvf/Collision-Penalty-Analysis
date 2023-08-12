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
// May 26, 2021 Theodore Kim, kim@cs.yale.edu
///////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef VOLUME_HYPERELASTIC_H
#define VOLUME_HYPERELASTIC_H

#include "SETTINGS.h"

namespace HOBAK {
namespace VOLUME {

class HYPERELASTIC
{
public:
  virtual ~HYPERELASTIC() = 0;

  // Computes the strain energy density
  virtual REAL psi(const MATRIX3& F) const;
  virtual REAL psi(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const;

  // Computes the first Piola-Kirchoff PK1 stress
  virtual MATRIX3 PK1(const MATRIX3& F) const;
  virtual MATRIX3 PK1(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const;

  // Computes the derivative of the PK1 stress
  virtual MATRIX9 hessian(const MATRIX3& F) const;
  virtual MATRIX9 hessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const; 

  // Computes the derivative of the PK1 stress, clamped to semi-positive definiteness
  virtual MATRIX9 clampedHessian(const MATRIX3& F) const;
  virtual MATRIX9 clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const;

  // The name of the material
  virtual std::string name() const = 0;

  // True if the energy computation requires the SVD of F
  virtual bool energyNeedsSVD() const = 0;

  // True if the PK1 computation requires the SVD of F
  virtual bool PK1NeedsSVD() const = 0;

  // convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
  static REAL computeMu(const REAL E, const REAL nu);
  static REAL computeLambda(const REAL E, const REAL nu);
};

}
}

#endif
