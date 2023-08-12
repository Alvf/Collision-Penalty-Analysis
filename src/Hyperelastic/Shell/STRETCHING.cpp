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

#include "STRETCHING.h"
#include "util/MATRIX_UTIL.h"

namespace HOBAK {
namespace SHELL {

STRETCHING::~STRETCHING() = default;

STRETCHING::STRETCHING(const REAL& mu, const REAL& lambda) :
  _mu(mu), _lambda(lambda)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Do a brute-force clamping
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 STRETCHING::clampedHessian(const MATRIX3x2& F) const
{
  const MATRIX6 H = hessian(F);
  return clampEigenvalues(H);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL STRETCHING::computeMu(const REAL E, const REAL nu)
{
  return E / (2.0 * (1.0 + nu));
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// convert Young's modulus (E) and Poisson's ratio (nu) to Lam\'{e} parameters
//
// Tamstorf claims that Lambda computation should be slightly different from the volume case in
// his dissertation, "Large scale simulation of cloth and hair with contact", 2016
//
// Here we use the formula from page 34, Equation 4.2. The original volumetric conversion
// is commented out.
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL STRETCHING::computeLambda(const REAL E, const REAL nu)
{
  //return (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));

  // 2016 version
  return (E * nu) / (1.0 - nu * nu);
}

} // SHELL
} // HOBAK
