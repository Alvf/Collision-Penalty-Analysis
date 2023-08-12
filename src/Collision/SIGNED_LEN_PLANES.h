/*
This is Alvin's attempt to streamline all of the chain rule stuff required
in calculating out penalty collision energies.
This object stores derivatives and values for plane-plane length
given some 3-edge configuration c.
*/

#ifndef SIGNED_LEN_PLANES_H
#define SIGNED_LEN_PLANES_H

#include "SETTINGS.h"
using namespace std;

class SIGNED_LEN_PLANES
{
public:
    SIGNED_LEN_PLANES(){}
    virtual ~SIGNED_LEN_PLANES() {};

    virtual string name()
    {
        return string("Signed Length (plane/plane)");
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

        return e2p.dot(e0p.cross(e1p))/e0p.cross(e1p).norm(); 
    }

    virtual VECTOR9 DlDc(VECTOR9 c) const
    {
        VECTOR3 e0p, e1p, e2p;
        getVecs(c, e0p, e1p, e2p);

        VECTOR3 cr01 = e0p.cross(e1p);
        VECTOR3 cr12 = e1p.cross(e2p);
        VECTOR3 cr20 = e2p.cross(e0p);

        REAL cr01n3 = cr01.norm()*cr01.norm()*cr01.norm();

        VECTOR3 DlDe0 = (e1p*e0p.dot(e1p)/e1p.dot(e1p) - e0p)*(e0p.dot(cr12))*(e1p.dot(e1p))/cr01n3 + cr12/cr01.norm();
        VECTOR3 DlDe1 = (e0p*e1p.dot(e0p)/e0p.dot(e0p) - e1p)*(e0p.dot(cr12))*(e0p.dot(e0p))/cr01n3 + cr20/cr01.norm();
        VECTOR3 DlDe2 = cr01/cr01.norm();

        VECTOR9 DlDc;
        DlDc.block<3,1>(0,0) = DlDe0;
        DlDc.block<3,1>(3,0) = DlDe1;
        DlDc.block<3,1>(6,0) = DlDe2;

        return DlDc;
    }

    virtual MATRIX9 d2lDc2(VECTOR9 c) const
    {
        VECTOR3 e0p, e1p, e2p;
        getVecs(c, e0p, e1p, e2p);
        // Eigenconstruction
        const REAL e0pN2 = e0p.squaredNorm();
        const REAL e1pN2 = e1p.squaredNorm();
        const REAL e2pN2 = e2p.squaredNorm();
        const REAL l = e2p.dot(e0p.cross(e1p))/e0p.cross(e1p).norm();
        
        const REAL f12 = sqrt(1 + 4*(e1pN2/e2pN2));
        const REAL f02 = sqrt(1 + 4*(e0pN2/e2pN2));

        const REAL laml0 = -l/(2*e1pN2)*(1 + f12);
        const REAL laml1 = -l/(2*e1pN2)*(1 - f12);
        const REAL laml2 = -l/(2*e0pN2)*(1 + f02);
        const REAL laml3 = -l/(2*e0pN2)*(1 - f02);

        MATRIX9 hess = MATRIX9::Zero();
        
        VECTOR9 q3 = VECTOR9::Zero();
        const REAL om3 = laml3/(laml3 - l/e2pN2);
        q3.block<3,1>(0,0) = e2p;
        q3.block<3,1>(6,0) = om3*e0p;
        q3.normalize();
        hess += laml3*q3*q3.transpose();

        VECTOR9 q2 = VECTOR9::Zero();
        const REAL om2 = laml2/(laml2 - l/e2pN2);
        q2.block<3,1>(0,0) = e2p;
        q2.block<3,1>(6,0) = om2*e0p;
        q2.normalize(); 
        hess += laml2*q2*q2.transpose();

        VECTOR9 q1 = VECTOR9::Zero();
        const REAL om1 = laml1/(laml1 - l/e2pN2);
        q1.block<3,1>(3,0) = e2p;
        q1.block<3,1>(6,0) = om1*e1p;
        q1.normalize();
        hess += laml1*q1*q1.transpose();

        VECTOR9 q0 = VECTOR9::Zero();
        const REAL om0 = laml0/(laml0 - l/e2pN2);
        q0.block<3,1>(3,0) = e2p;
        q0.block<3,1>(6,0) = om0*e1p;
        q0.normalize();
        hess += laml0*q0*q0.transpose();

        return hess;
    }

    virtual MATRIX9 clampedD2lDc(VECTOR9 c, REAL s = 1) const
    {
        VECTOR3 e0p, e1p, e2p;
        getVecs(c, e0p, e1p, e2p);
        // Eigenconstruction
        const REAL e0pN2 = e0p.squaredNorm();
        const REAL e1pN2 = e1p.squaredNorm();
        const REAL e2pN2 = e2p.squaredNorm();
        const REAL l = e2p.dot(e0p.cross(e1p))/e0p.cross(e1p).norm();
        
        const REAL f12 = sqrt(1 + 4*(e1pN2/e2pN2));
        const REAL f02 = sqrt(1 + 4*(e0pN2/e2pN2));

        const REAL laml0 = -l/(2*e1pN2)*(1 + f12);
        const REAL laml1 = -l/(2*e1pN2)*(1 - f12);
        const REAL laml2 = -l/(2*e0pN2)*(1 + f02);
        const REAL laml3 = -l/(2*e0pN2)*(1 - f02);

        MATRIX9 hess = MATRIX9::Zero();
        
        if(s*laml3 > 0.0){
            VECTOR9 q3 = VECTOR9::Zero();
            const REAL om3 = laml3/(laml3 - l/e2pN2);
            q3.block<3,1>(0,0) = e2p;
            q3.block<3,1>(6,0) = om3*e0p;
            q3.normalize();
            hess += s*laml3*q3*q3.transpose();
        }
        else{ //laml2 is the positive one...
            VECTOR9 q2 = VECTOR9::Zero();
            const REAL om2 = laml2/(laml2 - l/e2pN2);
            q2.block<3,1>(0,0) = e2p;
            q2.block<3,1>(6,0) = om2*e0p;
            q2.normalize(); 
            hess += s*laml2*q2*q2.transpose();
        }
        //lam1 or lam0?
        if(s*laml1 > 0.0){
            VECTOR9 q1 = VECTOR9::Zero();
            const REAL om1 = laml1/(laml1 - l/e2pN2);
            q1.block<3,1>(3,0) = e2p;
            q1.block<3,1>(6,0) = om1*e1p;
            q1.normalize();
            hess += s*laml1*q1*q1.transpose();
        }
        else{ //laml0 is the positive one...
            VECTOR9 q0 = VECTOR9::Zero();
            const REAL om0 = laml0/(laml0 - l/e2pN2);
            q0.block<3,1>(3,0) = e2p;
            q0.block<3,1>(6,0) = om0*e1p;
            q0.normalize();
            hess += s*laml0*q0*q0.transpose();
        }
        return hess;
    }

};

#endif
