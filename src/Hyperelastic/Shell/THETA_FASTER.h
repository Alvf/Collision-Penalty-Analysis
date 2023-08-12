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

#ifndef SHELL_THETA_FASTER_H
#define SHELL_THETA_FASTER_H

#include "THETA.h"

namespace HOBAK {
namespace SHELL {

class THETA_FASTER : public THETA
{
public:
  THETA_FASTER(const REAL& mu);
  virtual ~THETA_FASTER() { };

  // position-based Hessian
  virtual MATRIX12 hessian(const std::vector<VECTOR3>& flap, const REAL& restTheta) const;

  // The name of the material
  virtual std::string name() const;

protected:
  virtual VECTOR3 crossHessian0(const int i, const int j) const;
  virtual VECTOR3 crossHessian1(const int i, const int j) const;

  // populate _crossHessian0 and _crossHessian1
  void initializeCrossHessians();

  // cache out the cross Hessians
  VECTOR3 _crossHessian0[12][12];
  VECTOR3 _crossHessian1[12][12];
};

}
}

#endif
