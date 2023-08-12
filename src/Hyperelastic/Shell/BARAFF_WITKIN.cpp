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
//
// It is the shearing energy from the original Baraff-Witkin paper,
// "Large Steps in Cloth Simulation" SIGGRAPH 1998
//
// But, it's the FEM-style formulation and eigenanalysis from the Kim paper
// "A Finite Element Formulation of Baraff-Witkin Cloth" SCA 2020
///////////////////////////////////////////////////////////////////////////////////////////////////

#include "BARAFF_WITKIN.h"
#include "util/MATRIX_UTIL.h"

namespace HOBAK {
namespace SHELL {

BARAFF_WITKIN::BARAFF_WITKIN(const REAL& mu, const REAL& lambda) :
  STRETCHING(mu, lambda), _stretch(mu, lambda), _shear(mu, lambda)
{
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
std::string BARAFF_WITKIN::name() const
{ 
  return std::string("BARAFF_WITKIN"); 
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
REAL BARAFF_WITKIN::psi(const MATRIX3x2& F) const
{
  return _stretch.psi(F) + _shear.psi(F);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX3x2 BARAFF_WITKIN::PK1(const MATRIX3x2& F) const
{
  return _stretch.PK1(F) + _shear.PK1(F);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 BARAFF_WITKIN::hessian(const MATRIX3x2& F) const
{
  return _stretch.hessian(F) + _shear.hessian(F);
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
MATRIX6 BARAFF_WITKIN::clampedHessian(const MATRIX3x2& F) const
{
  // this clamp is conservative, so it won't match the brute force exactly
  return _stretch.clampedHessian(F) + _shear.clampedHessian(F);
}

} // SHELL
} // HOBAK
