/*
This is Alvin's effort to consolidate collision energies into 1D functions on general notions of length
This object is just a big parent.
*/

#ifndef ENERGY_1D_H
#define ENERGY_1D_H

#include "SETTINGS.h"

using namespace std;
//To save some space, energy1D itself will just be the 1D signed length spring function.
class ENERGY_1D
{
public:
  /*
    // some energies have different amounts of parameters. Just shove in those parameters as a vector and use them later.
    // TODO: looks like this only has two params -- should just make it explicit 
    ENERGY_1D(const vector<REAL> params)
    {
        _mu = params[0];
        _eps = params[1];
    };
    */
    
    ENERGY_1D(const REAL& mu, const REAL& collisionEps) :
      _mu(mu), _eps(collisionEps)
    {
    }

    virtual ~ENERGY_1D() {};
    const REAL& mu() const  { return _mu; };
    const REAL& eps() const { return _eps; };

    // redundant, and subject to corruption
    //vector<REAL> parameters;

    // just a nice descriptive name
    virtual string name(){
        return string("Signed Length Quadratic Penalty");
    }

    // fed some general length and some extra info at runtime, return the energy
    virtual REAL psi(const REAL& l, const vector<int>& extraInfo) const{
        int sgn = 1;
        if (extraInfo.size() > 0){
            sgn = extraInfo[0];
        }
        return _mu*(sgn*l - _eps)*(sgn*l - _eps);
    }
    virtual REAL psi(const REAL& l) const{
        vector<int> emptyExtras;
        return psi(l, emptyExtras);
    }

    // fed some general length and some extra info at runtime, return dpsi/dl
    virtual REAL DpsiDl(const REAL& l, const vector<int>& extraInfo) const{
        int sgn = 1;
        if (extraInfo.size() > 0){
            sgn = extraInfo[0];
        }
        return 2*_mu*(l - sgn*_eps);
    }
    virtual REAL DpsiDl(const REAL& l) const{
        vector<int> emptyExtras;
        return DpsiDl(l, emptyExtras);
    }

    // fed some general length and some extra info at runtime, return dpsi/dl
    virtual REAL D2psiDl2(const REAL& l, const vector<int>& extraInfo) const{
        return 2*_mu;
    }
    virtual REAL D2psiDl2(const REAL& l) const{
        vector<int> emptyExtras;
        return D2psiDl2(l, emptyExtras);
    }
protected:
    // stiffness
    REAL _mu;

    // collision epsilon
    REAL _eps;

};

#endif
