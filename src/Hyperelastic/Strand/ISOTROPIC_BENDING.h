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

#ifndef STRAND_ISOTROPIC_BENDING_H
#define STRAND_ISOTROPIC_BENDING_H

#include "SETTINGS.h"
#include "MATRIX_UTIL.h"

#include <iostream>

namespace HOBAK {
namespace STRAND {

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
class ISOTROPIC_BENDING
{
public:
  ISOTROPIC_BENDING(const REAL& mu) { _mu = mu; };
  ISOTROPIC_BENDING(const REAL& mu, const REAL& theta0) { _mu = mu; _theta0 = theta0; };
  virtual ~ISOTROPIC_BENDING() { };

  // Computes the strain energy density
  // E is the matrix of edge vectors, stacked column-wise
  virtual REAL psi(const MATRIX3x2& E, const REAL& theta0) const = 0;
  virtual MATRIX3x2 PK1(const MATRIX3x2& E, const REAL& theta0) const = 0;
  virtual MATRIX6 hessian(const MATRIX3x2& E, const REAL& theta0) const = 0;
  virtual MATRIX6 clampedHessian(const MATRIX3x2& E, const REAL& theta0) const = 0;
  virtual REAL psi(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const = 0;
  virtual MATRIX3x2 PK1(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const = 0;
  virtual MATRIX6 hessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const = 0;
  virtual MATRIX6 clampedHessian(const MATRIX3x2& E, const REAL& theta0, const bool inverted) const = 0;
  
  REAL psi(const MATRIX3x2& E) const                { return psi(E, _theta0); };
  MATRIX3x2 PK1(const MATRIX3x2& E) const           { return PK1(E, _theta0); };
  MATRIX6 hessian(const MATRIX3x2& E) const         { return hessian(E, _theta0); };
  MATRIX6 clampedHessian(const MATRIX3x2& E) const  { return clampedHessian(E, _theta0); };
  
  REAL psi(const MATRIX3x2& E, const bool inverted) const                { return psi(E, _theta0, inverted); };
  MATRIX3x2 PK1(const MATRIX3x2& E, const bool inverted) const           { return PK1(E, _theta0, inverted); };
  MATRIX6 hessian(const MATRIX3x2& E, const bool inverted) const         { return hessian(E, _theta0, inverted); };
  MATRIX6 clampedHessian(const MATRIX3x2& E, const bool inverted) const  { return clampedHessian(E, _theta0, inverted); };

  REAL& mu() { return _mu; };
  const REAL& mu() const { return _mu; };
  const REAL& theta0() const { return _theta0; };

  // The name of the material
  virtual std::string name() const = 0;

  // compute the angle between two edges formed by three particles
  static REAL angle(const VECTOR3& v0, const VECTOR3& v1, const VECTOR3& v2)
  {
    MATRIX3x2 E;
    E.col(0) = v0 - v1;
    E.col(1) = v2 - v1;
    return angle(E);
  }

  // compute the angle between two edges
  static REAL angle(const MATRIX3x2& E)
  {
    const VECTOR3 e0 = E.col(0);
    const VECTOR3 e1 = E.col(1);
    //const REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
    REAL cosTheta = e0.dot(e1) / (e0.norm() * e1.norm());
    cosTheta = (cosTheta > 1.0) ? 1.0 : cosTheta;
    cosTheta = (cosTheta < -1.0) ? -1.0 : cosTheta;
    return acos(cosTheta);
  }

protected:
  REAL _mu;
  REAL _theta0;
};

}
}

#endif
