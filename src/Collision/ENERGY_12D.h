/*
psi(l(c(x))) and dpsi/dx with chain rule stuff and interchangable l(c)s, c(x)s, and psi(l)s
*/

#ifndef ENERGY_12D_H
#define ENERGY_12D_H

#include <iostream>

#include "SETTINGS.h"

#include "ENERGY_1D.h"
#include "C_PLANES.h"
#include "SIGNED_LEN_PLANES.h"

// DEBUG
#include "MATRIX_UTIL.h"

using namespace std;

class ENERGY_12D
{
public:
    ENERGY_12D(ENERGY_1D* psi, C_PLANES* c, SIGNED_LEN_PLANES* l) :
    psiFunc(psi), frameFunc(c), lenFunc(l) {}
    virtual ~ENERGY_12D() { };

    ENERGY_1D* psiFunc;
    C_PLANES* frameFunc;
    SIGNED_LEN_PLANES* lenFunc;

    // stiffness
    const REAL& mu() const { return psiFunc->mu(); };

    // collision epsilon
    const REAL& eps() const { return psiFunc->eps(); };

    virtual string name() const {
        string psiName = psiFunc->name();
        string frameName = frameFunc->name();
        string lenName = lenFunc->name();
        return psiName + string(" | ") + frameName + string(" | ") + lenName;
    }

    virtual VECTOR12 makeState(const vector<VECTOR3>& vs) const
    {
        VECTOR12 x = VECTOR12::Zero();
        for(int i = 0; i < 4; i++){
            x.block<3,1>(3*i,0) = vs[i];
        }
        return x;
    }

    virtual REAL getLen(const VECTOR12& x) const
    {
        return lenFunc->len(frameFunc->getC(x));
    }

    virtual REAL psi(const VECTOR12& x, const vector<int>& extraInfo) const
    {
        const REAL l = lenFunc->len(frameFunc->getC(x));
        if(isnan(l)) return 0;
        return psiFunc->psi(l, extraInfo);
    }
    REAL psi(const VECTOR12& x) const
    {
        vector<int> emptyArgs;
        return psi(x, emptyArgs);
    }
    REAL psi(const vector<VECTOR3>& vs) const
    {
        vector<int> emptyArgs;
        return psi(makeState(vs),emptyArgs);
    }
    REAL psi(const vector<VECTOR3>& vs, const vector<int>& extraInfo) const
    {
        return psi(makeState(vs), extraInfo);
    }

    virtual VECTOR12 gradient(const VECTOR12& x, const vector<int>& extraInfo) const
    {
        const VECTOR9 c = frameFunc->getC(x);
        const REAL l = lenFunc->len(c);
        if (isnan(l)) return VECTOR12::Zero();

        const REAL dpdl = psiFunc->DpsiDl(l, extraInfo);
        const VECTOR9 gl = lenFunc->DlDc(c);
        const MATRIX9x12 dcdx = frameFunc->dCdx_s(x);

        return dpdl*dcdx.transpose()*gl;
    }
    VECTOR12 gradient(const VECTOR12& x) const
    {
        vector<int> emptyArgs;
        return gradient(x, emptyArgs);
    }
    VECTOR12 gradient(const vector<VECTOR3>& vs) const
    {
        vector<int> emptyArgs;
        return gradient(makeState(vs), emptyArgs);
    }
    VECTOR12 gradient(const vector<VECTOR3>& vs, const vector<int>& extraInfo) const
    {
        return gradient(makeState(vs), extraInfo);
    }

    virtual MATRIX12 hessian(const VECTOR12& x, const vector<int>& extraInfo) const
    {
        const VECTOR9 c = frameFunc->getC(x);
        const REAL l = lenFunc->len(c);
        if(isnan(l)) return MATRIX12::Zero();

        const REAL dpdl = psiFunc->DpsiDl(l, extraInfo);
        const REAL d2pdl2 = psiFunc->D2psiDl2(l, extraInfo);
        const VECTOR9 gl = lenFunc->DlDc(c);
        const MATRIX9 Hl = lenFunc->d2lDc2(c);
        const MATRIX9x12 dcdx = frameFunc->dCdx_s(x);
        
        if (!frameFunc->has_del()){
            return dcdx.transpose() * ( d2pdl2*gl*gl.transpose() + dpdl*Hl ) * dcdx;
        }
        const MATRIX9x12 del_dcdx = frameFunc->del_dCdx_s(x);
        return dcdx.transpose() * (d2pdl2*gl*gl.transpose() + dpdl*Hl) * dcdx - del_dcdx.transpose()*dpdl*Hl*del_dcdx;
    }
    MATRIX12 hessian(const VECTOR12& x) const
    {
        vector<int> emptyArgs;
        return hessian(x, emptyArgs);
    }
    MATRIX12 hessian(const vector<VECTOR3>& vs) const
    {
        vector<int> emptyArgs;
        return hessian(makeState(vs), emptyArgs);
    }
    MATRIX12 hessian(const vector<VECTOR3>& vs, const vector<int>& extraInfo) const
    {
        return hessian(makeState(vs), extraInfo);
    }

    virtual MATRIX12 clampedHessian(const VECTOR12& x, const vector<int>& extraInfo) const
    {
        const VECTOR9 c = frameFunc->getC(x);
        const REAL l = lenFunc->len(c);
        if(isnan(l)) return MATRIX12::Zero();

        const REAL dpdl = psiFunc->DpsiDl(l, extraInfo);
        const REAL d2pdl2 = psiFunc->D2psiDl2(l, extraInfo);

        MATRIX9 clampedD2Pdc2 = MATRIX9::Zero();

        const MATRIX9 clampedHl = lenFunc->clampedD2lDc(c,dpdl);

        clampedD2Pdc2 += clampedHl;

        if(d2pdl2>0.0)
        {
            VECTOR9 gl = lenFunc->DlDc(c);
            gl.normalize();
            clampedD2Pdc2 += d2pdl2*gl*gl.transpose();
        }

        const MATRIX9x12 dcdx = frameFunc->dCdx_s(x);

        if(!frameFunc->has_del()){
            return dcdx.transpose()*clampedD2Pdc2*dcdx;
        }
        MATRIX9x12 del_dcdx = frameFunc->del_dCdx_s(x);
        MATRIX9 negClamped = lenFunc->clampedD2lDc(c, -dpdl);
        return dcdx.transpose()*clampedD2Pdc2*dcdx + del_dcdx.transpose()*negClamped*del_dcdx;
    }
    MATRIX12 clampedHessian(const VECTOR12& x) const
    {
        vector<int> emptyArgs;
        return clampedHessian(x, emptyArgs);
    }
    MATRIX12 clampedHessian(const vector<VECTOR3>& vs) const
    {
        vector<int> emptyArgs;
        return clampedHessian(makeState(vs), emptyArgs);
    }
    MATRIX12 clampedHessian(const vector<VECTOR3>& vs, const vector<int>& extraInfo) const
    {
        return clampedHessian(makeState(vs), extraInfo);
    }
};

#endif
