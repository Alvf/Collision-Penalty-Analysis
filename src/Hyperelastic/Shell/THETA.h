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

#ifndef SHELL_THETA_H
#define SHELL_THETA_H

#include "BENDING.h"

namespace HOBAK {
namespace SHELL {

class THETA : public BENDING
{
public:
  THETA(const REAL& mu);
  virtual ~THETA() { };

  // Computes the strain energy density
  virtual REAL psi(const std::vector<VECTOR3>& flap, const REAL& restTheta) const;

  // position-based gradient
  virtual VECTOR12 gradient(const std::vector<VECTOR3>& flap, const REAL& restTheta) const;

  // position-based Hessian
  virtual MATRIX12 hessian(const std::vector<VECTOR3>& flap, const REAL& restTheta) const;

  // The name of the material
  virtual std::string name() const;

  // Compute the bending angle based on a rest flap.
  virtual REAL restAngle(const std::vector<VECTOR3>& restFlap) const;

protected:
  VECTOR3 normal0(const std::vector<VECTOR3>& flap) const;
  VECTOR3 normal1(const std::vector<VECTOR3>& flap) const;
  MATRIX3x12 normalGradient0(const std::vector<VECTOR3>& flap) const;
  MATRIX3x12 normalGradient1(const std::vector<VECTOR3>& flap) const;
  virtual std::vector<MATRIX12> normalHessian0(const std::vector<VECTOR3>& flap) const;
  virtual std::vector<MATRIX12> normalHessian1(const std::vector<VECTOR3>& flap) const;
  
  MATRIX3x12 edgeGradient(const std::vector<VECTOR3>& flap) const;

  REAL cosThetaPsi(const std::vector<VECTOR3>& flap) const;
  VECTOR12 cosThetaGradient(const std::vector<VECTOR3>& flap) const;
  MATRIX12 cosThetaHessian(const std::vector<VECTOR3>& flap) const;
  
  REAL sinThetaPsi(const std::vector<VECTOR3>& flap) const;
  VECTOR12 sinThetaGradient(const std::vector<VECTOR3>& flap) const;
  MATRIX12 sinThetaHessian(const std::vector<VECTOR3>& flap) const;

  virtual VECTOR3 crossHessian0(const int i, const int j) const;
  virtual VECTOR3 crossHessian1(const int i, const int j) const;

  static MATRIX3x12 crossGradient0(const std::vector<VECTOR3>& flap);
  static MATRIX3x12 crossGradient1(const std::vector<VECTOR3>& flap);
};

}
}

#endif
