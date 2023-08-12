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
#ifndef VOLUME_ARAP_DAMPING_H
#define VOLUME_ARAP_DAMPING_H

#include "DAMPING.h"

namespace HOBAK {
namespace VOLUME {

class ARAP_DAMPING : public DAMPING
{
public:
    ARAP_DAMPING(const REAL& mu);
    ~ARAP_DAMPING() {};

    virtual REAL psi(const MATRIX3& F, const MATRIX3& Fdot) const;

    virtual MATRIX3 PK1(const MATRIX3& F, const MATRIX3& Fdot) const;

    virtual std::string name() const;

    virtual MATRIX9 hessian(const MATRIX3& F, const MATRIX3& Fdot) const;

    virtual MATRIX9 clampedHessian(const MATRIX3& F, const MATRIX3& Fdot) const;
    
};

}
}

#endif
