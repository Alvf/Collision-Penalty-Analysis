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
#include "GREEN_DAMPING.h"
#include "MATRIX_UTIL.h"

namespace HOBAK {
namespace VOLUME {

GREEN_DAMPING::GREEN_DAMPING(const REAL& mu)
{
  _mu = mu;
}

std::string GREEN_DAMPING::name() const
{
  return "GREEN_DAMPING";
}

///////////////////////////////////////////////////////////////////////
// get the strain energy
///////////////////////////////////////////////////////////////////////
REAL GREEN_DAMPING::psi(const MATRIX3& F, const MATRIX3& Fdot) const
{
  const MATRIX3 FdotF = Fdot.transpose() * F;
  const MATRIX3 Edot = 0.5 * (FdotF + FdotF.transpose());
  return _mu * Edot.squaredNorm();
}

///////////////////////////////////////////////////////////////////////
// P = first Piola-Kirchoff stress tensor
///////////////////////////////////////////////////////////////////////
MATRIX3 GREEN_DAMPING::PK1(const MATRIX3& F, const MATRIX3& Fdot) const
{
  return _mu * F * (F.transpose() * Fdot + Fdot.transpose() * F);
}

///////////////////////////////////////////////////////////////////////
// there's probably a cleaner outer product form to be discovered
// here
///////////////////////////////////////////////////////////////////////
MATRIX9 GREEN_DAMPING::hessian(const MATRIX3& F, const MATRIX3& Fdot) const
{
  MATRIX9 pPpF = MATRIX9::Zero();
  int index = 0; 
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++, index++)
    {
      MATRIX3 pFpF;
      partialFpartialF(i, j, pFpF);
      MATRIX3 column =  F * F.transpose() * pFpF + F * pFpF.transpose() * F;
      pPpF.col(index) = flatten(column);
    }

  return _mu * pPpF;
}

///////////////////////////////////////////////////////////////////////
// filtered energy hessian
///////////////////////////////////////////////////////////////////////
MATRIX9 GREEN_DAMPING::clampedHessian(const MATRIX3& F, const MATRIX3& Fdot) const
{
  // analysis shows that it can never be indefinite, so no need to
  // filter
  //return hessian(F, Fdot);

  // doing some debugging ...
  return clampEigenvalues(hessian(F, Fdot));
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX9 GREEN_DAMPING::positionGradient(const MATRIX3& F, const MATRIX3& Fdot) const
{
  MATRIX9 pPpF = MATRIX9::Zero();
  int index = 0;
  MATRIX3 Fsum = F.transpose() * Fdot + Fdot.transpose() * F; 
  MATRIX3 Fproduct = F * Fdot.transpose();
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++, index++)
    {
      MATRIX3 pFpF;
      partialFpartialF(i, j, pFpF);
      MATRIX3 column =  pFpF * Fsum + F * pFpF.transpose() * Fdot + Fproduct * pFpF;
      pPpF.col(index) = flatten(column);
    }

  return _mu * pPpF;
}

}
}
