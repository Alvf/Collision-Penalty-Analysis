/*
This is Alvin's effort to consolidate collision energies into 1D functions on general notions of length
This object is just a big parent.
*/

#ifndef CUBIC_SIGNED_SPRING_1D_H
#define CUBIC_SIGNED_SPRING_1D_H

#include "ENERGY_1D.h"

using namespace std;
//To save some space, energy1D itself will just be the 1D signed length spring function.
class CUBIC_SIGNED_SPRING_1D : public ENERGY_1D
{
public:
    /* 
    // some energies have different amounts of parameters. Just shove in those parameters as a vector and use them later.
    CUBIC_SIGNED_SPRING1D(const vector<REAL> params) : ENERGY_1D(params)
    {
        _del = params[2];
    };
    */

    CUBIC_SIGNED_SPRING_1D(const REAL& mu, const REAL& collisionEps, const REAL& del) : 
      ENERGY_1D(mu, collisionEps), _del(del)
    {
    }

    CUBIC_SIGNED_SPRING_1D(const REAL& mu, const REAL& collisionEps) : 
      ENERGY_1D(mu, collisionEps)
    {
        _del = 0.001;
    }

    // just a nice descriptive name
    virtual string name() override{
        return string("Signed Length Cubic Penalty");
    }

    // fed some general length, return the energy
    virtual REAL psi(const REAL& l, const vector<int>& extraInfo) const override{
        REAL sQsTerm = sqrt((l-_eps)*(l-_eps) + _del);
        return _mu*(l-_eps)*(l-_eps)*sQsTerm;

    }

    // fed some general length, return dpsi/dl
    virtual REAL DpsiDl(const REAL& l, const vector<int>& extraInfo) const override{
        REAL sQsTerm = sqrt((l-_eps)*(l-_eps) + _del);
        return _mu*(l-_eps)*(2*_del + 3*(l-_eps)*(l-_eps))/sQsTerm;
    }

    // fed some general length, return dpsi/dl
    virtual REAL D2psiDl2(const REAL& l, const vector<int>& extraInfo) const override{
        REAL OME2 = (l-_eps)*(l-_eps);
        REAL sTerm = OME2 + _del;
        return _mu*(2*_del*_del + 9*_del*OME2 + 6*OME2*OME2)/(sTerm*sqrt(sTerm));
    }

private:
    REAL _del;
};

#endif
