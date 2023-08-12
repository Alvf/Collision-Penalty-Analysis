
/*
This is Alvin's attempt to streamline all of the chain rule stuff required
in calculating out penalty collision energies.
This object stores derivatives and values for a special C_ee, the frame used in the plane-plane case
given some 4-point configuration x (R12).
WARNING: A1 + A2 DNE -A1 when using this frame! Needs numerical eigenclamping, or just 
be fine with the Hessian being wrong (seems to look fine)
*/

#ifndef C_SPECIAL_EE_H
#define C_SPECIAL_EE_H

#include "SETTINGS.h"
#include "C_PLANES.h"
#include "Geometry/LINE_INTERSECT.h"
using namespace std;

class C_SPECIAL_EE : public C_PLANES
{
    public:
    C_SPECIAL_EE() : C_PLANES() {}
    virtual ~C_SPECIAL_EE() {};

    virtual string name() const override
    {
        return string("Collision Frame for special unsigned edge-edge collisions");
    }
    
    virtual VECTOR12 dal_dx(const VECTOR12& x) const override{
        vector<VECTOR3> v;
        v.resize(4);
        for (int i = 0; i < 4; i++)
        {
            v[i][0] = x[i * 3];
            v[i][1] = x[i * 3 + 1];
            v[i][2] = x[i * 3 + 2];
        }

        VECTOR3 e0, e1, e2;

        e0 = v[1] - v[0];
        e1 = v[3] - v[2];
        e2 = v[2] - v[0];

        const REAL e0d0 = e0.dot(e0);
        const REAL e0d1 = e0.dot(e1);
        const REAL e0d2 = e0.dot(e2);
        const REAL e1d1 = e1.dot(e1);
        const REAL e1d2 = e1.dot(e2);
        const REAL e0x12 = (e0.cross(e1)).squaredNorm();

        VECTOR12 dadx = VECTOR12::Zero();
        REAL al = (e0d2*e1d1 - e1d2 * e0d1)/(e0x12);
        if (al < 0 || al > 1) return dadx;
        
        const VECTOR3 dadv0 = e0d1 * e1 + e1d2 * e1 - e1d1 * (e0 + e2) - 2*al*(e1*e0d1 - e0*e1d1);
        const VECTOR3 dadv1 = e1d1*e2 - e1d2*e1 + 2*al*(e1*e0d1 - e0*e1d1);
        const VECTOR3 dadv2 = e1d2*e0 + e1d1*e0 - 2*e0d2*e1 - e0d1*(e1 - e2) - 2*al*(e0*e0d1 - e1*e0d0);
        const VECTOR3 dadv3 = 2*e0d2*e1 - e0d1*e2 - e1d2*e0 + 2*al*(e0*e0d1 - e1*e0d0);

        dadx.block<3,1>(0,0) = dadv0;
        dadx.block<3,1>(3,0) = dadv1;
        dadx.block<3,1>(6,0) = dadv2;
        dadx.block<3,1>(9,0) = dadv3;

        return dadx/e0x12;
    }
    virtual VECTOR12 dbet_dx(const VECTOR12& x) const override{
        vector<VECTOR3> v;
        v.resize(4);
        for (int i = 0; i < 4; i++)
        {
            v[i][0] = x[i * 3];
            v[i][1] = x[i * 3 + 1];
            v[i][2] = x[i * 3 + 2];
        }

        VECTOR3 e0, e1, e2;

        e0 = v[1] - v[0];
        e1 = v[3] - v[2];
        e2 = v[2] - v[0];

        const REAL e0d0 = e0.dot(e0);
        const REAL e0d1 = e0.dot(e1);
        const REAL e0d2 = e0.dot(e2);
        const REAL e1d1 = e1.dot(e1);
        const REAL e1d2 = e1.dot(e2);
        const REAL e0x12 = (e0.cross(e1)).squaredNorm();

        REAL bet = (e0d1 * e0d2 - e1d2*e0d0)/(e0x12);
        VECTOR12 dbdx = VECTOR12::Zero();
        if (bet < 0 || bet > 1) return dbdx;
        
        const VECTOR3 dbdv0 = 2*e1d2*e0 - e0d2*e1 + e0d0*e1 - e0d1*(e2 + e0) - 2*bet*(e1*e0d1 - e0*e1d1);
        const VECTOR3 dbdv1 = -2*e1d2*e0 + e0d2*e1 + e0d1*e2 + 2*bet*(e1*e0d1 - e0*e1d1);
        const VECTOR3 dbdv2 = e0d1*e0 - e0d2*e0 - e0d0*(e1-e2) - 2*bet*(e0*e0d1 - e1*e0d0);
        const VECTOR3 dbdv3 = e0d2*e0 - e0d0*e2 + 2*bet*(e0*e0d1 - e1*e0d0);

        dbdx.block<3,1>(0,0) = dbdv0;
        dbdx.block<3,1>(3,0) = dbdv1;
        dbdx.block<3,1>(6,0) = dbdv2;
        dbdx.block<3,1>(9,0) = dbdv3;

        return dbdx/e0x12;
    }

    virtual VECTOR9 getC(const VECTOR12& x) const //x fed in as in paper
    {
        VECTOR9 c;
        vector<VECTOR3> v;
        v.resize(4);
        for (int i = 0; i < 4; i++)
        {
            v[i][0] = x[i * 3];
            v[i][1] = x[i * 3 + 1];
            v[i][2] = x[i * 3 + 2];
        }

        VECTOR3 e0, e1, e2;

        e0 = v[1] - v[0];
        e1 = v[3] - v[2];
        e2 = v[2] - v[0]; 

        const REAL e0d1 = e0.dot(e1);
        const REAL e0d2 = e0.dot(e2);
        const REAL e1d2 = e1.dot(e2);
        const REAL e0d0 = e0.dot(e0);
        const REAL e1d1 = e1.dot(e1);
        const REAL e2d2 = e2.dot(e2);
        const REAL e0x12 = e0.cross(e1).squaredNorm();

        // VECTOR3 outer, inner;
        // IntersectLineSegments(v[0],v[1],v[2],v[3],outer, inner);
        // REAL a = (outer-v[0]).norm()/e0.norm();
        // REAL b = (inner - v[2]).norm()/e1.norm(); 
        REAL a = (e0d2*e1d1 - e1d2*e0d1)/(e0x12);
        if (a < 0) a = 0;
        if (a > 1) a = 1;
        REAL b = (e0d1*e0d2 - e1d2*e0d0)/(e0x12);
        if (b < 0) b = 0;
        if (b > 1) b = 1;

        c.block<3,1>(0,0) = e0;
        c.block<3,1>(3,0) = e1;
        c.block<3,1>(6,0) = e2 - a*e0 + b*e1;

        return c;
    }

    virtual MATRIX9x12 dCdx_s(const VECTOR12& x) const
    {
        vector<VECTOR3> v;
        v.resize(4);
        for (int i = 0; i < 4; i++)
        {
            v[i][0] = x[i * 3];
            v[i][1] = x[i * 3 + 1];
            v[i][2] = x[i * 3 + 2];
        }

        VECTOR3 e0, e1, e2;

        e0 = v[1] - v[0];
        e1 = v[3] - v[2];
        e2 = v[2] - v[0];

        const REAL e0d1 = e0.dot(e1);
        const REAL e0d2 = e0.dot(e2);
        const REAL e1d2 = e1.dot(e2);
        const REAL e0d0 = e0.dot(e0);
        const REAL e1d1 = e1.dot(e1);
        const REAL e2d2 = e2.dot(e2);
        const REAL e0x12 = e0.cross(e1).squaredNorm();

        VECTOR3 outer, inner;
        REAL a = (e0d2*e1d1 - e1d2*e0d1)/(e0x12);
        if (a < 0) a = 0;
        if (a > 1) a = 1;
        REAL b = (e0d1*e0d2 - e1d2*e0d0)/(e0x12);
        if (b < 0) b = 0;
        if (b > 1) b = 1;

        MATRIX9x12 partFrame;
        MATRIX3 id = MATRIX3::Identity();
        partFrame.setZero();
        partFrame.block<3,3>(0,0) = -id;
        partFrame.block<3,3>(0,3) = id;
        partFrame.block<3,3>(3,6) = -id;
        partFrame.block<3,3>(3,9) = id;
        partFrame.block<3,3>(6,0) = (a-1)*id;
        partFrame.block<3,3>(6,3) = -a*id;
        partFrame.block<3,3>(6,6) = (1-b)*id;
        partFrame.block<3,3>(6,9) = b*id;

        return partFrame;
    }

    // del term is nonzero
    virtual bool has_del() override{
        return true;
    }

    //del term a bit different from others
    virtual MATRIX9x12 del_dCdx_s(const VECTOR12& x) const override
    {
        vector<VECTOR3> v;
        v.resize(4);
        for (int i = 0; i < 4; i++)
        {
            v[i][0] = x[i * 3];
            v[i][1] = x[i * 3 + 1];
            v[i][2] = x[i * 3 + 2];
        }

        VECTOR3 e0, e1, e2;

        e0 = v[1] - v[0];
        e1 = v[3] - v[2];
        e2 = v[2] - v[0];

        VECTOR12 dadx, dbdx, dgdx;
        dadx = dal_dx(x);
        dbdx = dbet_dx(x);

        MATRIX9x12 partFrame;
        partFrame.setZero();
        partFrame.block<3,12>(6,0) = -e0*dadx.transpose() + e1*dbdx.transpose();

        return partFrame;
    }

};

#endif
