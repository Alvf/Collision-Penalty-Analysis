/*
This is Alvin's effort to consolidate collision energies into 1D functions on general notions of length
This object is just a big parent.
*/

#ifndef UNSIGNED_SPRING_1D_H
#define UNSIGNED_SPRING_1D_H

#include "ENERGY_1D.h"

using namespace std;
//To save some space, energy1D itself will just be the 1D signed length spring function.
class UNSIGNED_SPRING_1D : public ENERGY_1D
{
public:
    /* 
    // some energies have different amounts of parameters. Just shove in those parameters as a vector and use them later.
    UNSIGNED_SPRING_1D(const vector<REAL> params) : ENERGY_1D(params)
    {
    };
    */

    UNSIGNED_SPRING_1D(const REAL& mu, const REAL& collisionEps) : ENERGY_1D(mu, collisionEps)
    {
    }

    // just a nice descriptive name
    virtual string name() override{
        return string("Unsigned Length Quadratic Penalty");
    }

    // fed some general length, return the energy
    virtual REAL psi(const REAL& l, const vector<int>& extraInfo) const override{
        REAL sgn = 1;
        if (extraInfo.size() !=0){
            sgn = sgn*extraInfo[0];
        }
        return _mu*(sgn*l - _eps)*(sgn*l - _eps);
    }

    // fed some general length, return dpsi/dl
    virtual REAL DpsiDl(const REAL& l, const vector<int>& extraInfo) const override{
        // extraInfo could contain no elements
        //REAL sgn = extraInfo[0] == 1 ? -1: 1;
        REAL sgn = 1;
        if (extraInfo.size() != 0) 
          sgn = sgn*extraInfo[0];
        return 2*_mu*(l - sgn*_eps);
    }

    // fed some general length, return dpsi/dl
    virtual REAL D2psiDl2(const REAL& l, const vector<int>& extraInfo) const override{
        return 2*_mu;
    }
};

#endif
