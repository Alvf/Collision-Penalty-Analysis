/*
This is Alvin's attempt to streamline all of the chain rule stuff required
in calculating out penalty collision energies.
This object stores derivatives and values for plane-plane length
given some 3-edge configuration c.
*/

#ifndef UNSIGNED_LEN_PLANES_H
#define UNSIGNED_LEN_PLANES_H

#include "SETTINGS.h"
#include "SIGNED_LEN_PLANES.h"
using namespace std;

class UNSIGNED_LEN_PLANES : public SIGNED_LEN_PLANES
{
public:
    UNSIGNED_LEN_PLANES(): SIGNED_LEN_PLANES(){}
    virtual ~UNSIGNED_LEN_PLANES() {};

    virtual string name()
    {
        return string("Unsigned Length (plane/plane)");
    }

    virtual void getVecs(VECTOR9 c, VECTOR3& e0p, VECTOR3& e1p, VECTOR3& e2p) const
    {
        e0p = c.block<3,1>(0,0);
        e1p = c.block<3,1>(3,0);
        e2p = c.block<3,1>(6,0);
    }

    virtual REAL len(VECTOR9 c) const
    {
        VECTOR3 e0p, e1p, e2p;
        getVecs(c, e0p, e1p, e2p);

        return e2p.norm(); 
    }

    virtual VECTOR9 DlDc(VECTOR9 c) const
    {
        VECTOR3 e0p, e1p, e2p;
        getVecs(c, e0p, e1p, e2p);

        VECTOR9 DlDc = VECTOR9::Zero();
        DlDc.block<3,1>(6,0) = e2p.normalized();

        return DlDc;
    }

    virtual MATRIX9 d2lDc2(VECTOR9 c) const
    {
        VECTOR3 e0p, e1p, e2p;
        getVecs(c, e0p, e1p, e2p);
        // Eigenconstruction
        const REAL u = e2p.norm();

        MATRIX9 hess = MATRIX9::Zero();
        
        MATRIX3 projBlock = MATRIX3::Identity() - e2p*e2p.transpose()/(u*u);
        hess.block<3,3>(6,6) = projBlock/u;

        return hess;
    }

    virtual MATRIX9 clampedD2lDc(VECTOR9 c, REAL s = 1) const
    {
        VECTOR3 e0p, e1p, e2p;
        getVecs(c, e0p, e1p, e2p);
        // Eigenconstruction
        const REAL u = e2p.norm();
        
        const REAL laml0 = 1/u;

        MATRIX9 hess = MATRIX9::Zero();
        if (s*laml0 < 0.0){
            return hess;
        }
        
        MATRIX3 projBlock = MATRIX3::Identity() - e2p*e2p.transpose()/(u*u);
        hess.block<3,3>(6,6) = s*laml0*projBlock;

        return hess;
    }

};

#endif
