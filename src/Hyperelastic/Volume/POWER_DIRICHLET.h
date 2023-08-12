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
#ifndef VOLUME_POWER_DIRICHLET_H
#define VOLUME_POWER_DIRICHLET_H

#include "HYPERELASTIC.h"

namespace HOBAK {
namespace VOLUME {

class POWER_DIRICHLET : public HYPERELASTIC
{
public:
    POWER_DIRICHLET(const REAL& mu, const VECTOR3& a);
    ~POWER_DIRICHLET() {};

    // get the strain energy
    virtual REAL psi(const MATRIX3& F) const override;
    virtual MATRIX3 PK1(const MATRIX3& F) const override;

    virtual std::string name() const override;

    virtual MATRIX9 hessian(const MATRIX3& F) const override;

    virtual MATRIX9 clampedHessian(const MATRIX3& F) const override;
    
    virtual bool energyNeedsSVD() const override;

    virtual bool PK1NeedsSVD() const override;

private:
    REAL _mu;
    REAL _power;

    VECTOR3 _a;
    MATRIX3 _A;
};

} // VOLUME
} // HOBAK

#endif
