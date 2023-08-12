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
#ifndef VOLUME_SNH_H
#define VOLUME_SNH_H

#include "HYPERELASTIC.h"

namespace HOBAK {
namespace VOLUME {

class SNH final : public HYPERELASTIC
{
public:
  SNH(const REAL& mu, const REAL& lambda);

  virtual ~SNH() override = default;

  virtual REAL psi(const MATRIX3& F) const override;

  virtual MATRIX3 PK1(const MATRIX3& F) const override;

  virtual MATRIX9 hessian(const MATRIX3& F) const override;

  virtual MATRIX9 clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const override;

  virtual std::string name() const override;

  virtual bool energyNeedsSVD() const override;

  virtual bool PK1NeedsSVD() const override;

private:
  const REAL _mu;
  const REAL _lambda;
  const REAL _alpha;
};

} // VOLUME
} // HOBAK

#endif
