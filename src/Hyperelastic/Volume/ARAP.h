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
#ifndef VOLUME_ARAP_H
#define VOLUME_ARAP_H

#include "HYPERELASTIC.h"

namespace HOBAK {
namespace VOLUME {

class ARAP : public HYPERELASTIC
{
public:
    ARAP(const REAL& mu, const REAL& lambda);
    ~ARAP() {};

    // get the strain energy
    virtual REAL psi(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const override;

    virtual MATRIX3 PK1(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const override;

    virtual std::string name() const override;

    virtual MATRIX9 hessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const override;

    virtual MATRIX9 clampedHessian(const MATRIX3& U, const VECTOR3& Sigma, const MATRIX3& V) const override;
    
    virtual bool energyNeedsSVD() const override;

    virtual bool PK1NeedsSVD() const override;

private:
    REAL _lambda;
    REAL _mu;
};

} // VOLUME
} // HOBAK

#endif
