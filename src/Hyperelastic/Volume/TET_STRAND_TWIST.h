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
#ifndef VOLUME_TET_STRAND_TWIST_H
#define VOLUME_TET_STRAND_TWIST_H

#include "HYPERELASTIC.h"

namespace HOBAK {
namespace VOLUME {

class TET_STRAND_TWIST : public HYPERELASTIC
{
public:
    TET_STRAND_TWIST(const REAL& mu, const REAL& theta0);
    ~TET_STRAND_TWIST() {};

    // get the strain energy
    virtual REAL psi(const MATRIX3& F) const override;

    virtual MATRIX3 PK1(const MATRIX3& F) const override;

    virtual std::string name() const override;

    virtual MATRIX9 hessian(const MATRIX3& F) const override;
    
    virtual MATRIX9 clampedHessian(const MATRIX3& F) const override;
    
    MATRIX12 forceGradient(const std::vector<VECTOR3>& vertices) const;
    MATRIX12 clampedForceGradient(const std::vector<VECTOR3>& vertices) const;
    
    MATRIX12 forceGradientRank4(const std::vector<VECTOR3>& vertices) const;
    MATRIX12 clampedForceGradientRank4(const std::vector<VECTOR3>& vertices) const;

    virtual bool energyNeedsSVD() const override;

    virtual bool PK1NeedsSVD() const override;

    const REAL& theta0() const { return _theta0; };
    REAL& theta0()             { return _theta0; };

    REAL& mu()             { return _mu; };
    const REAL& mu() const { return _mu; };

private:
    static void computeAlphaBeta(const std::vector<VECTOR3>& vertices, REAL& alpha, REAL& beta);
    static MATRIX3 computeFtwist(const std::vector<VECTOR3>& vertices);
    static MATRIX9x12 computeTwistPFPx(const std::vector<VECTOR3>& vertices);
    static MATRIX9x12 computeDelTwistPFPx(const std::vector<VECTOR3>& vertices);
    static MATRIX9x12 computeTwistPFPxFull(const std::vector<VECTOR3>& vertices);
    MATRIX12 computeRank4Correction(const std::vector<VECTOR3>& vertices) const;
    MATRIX12 computeRank4CorrectionClamped(const std::vector<VECTOR3>& vertices) const;

    static void computeFandGradients(const std::vector<VECTOR3>& vertices, 
                                     MATRIX3& Ftwist, MATRIX9x12& pFpX, MATRIX9x12& delpFpX);
    void clampedHessian(const MATRIX3& F, MATRIX9& plusH, MATRIX9& minusH) const;

    REAL _mu;
    REAL _theta0;
    REAL _zeroGuard;
};

} // VOLUME
} // HOBAK

#endif
