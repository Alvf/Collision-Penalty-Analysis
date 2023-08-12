
/*
This is Alvin's attempt to streamline all of the chain rule stuff required
in calculating out penalty collision energies.
This object stores derivatives and values for C_vf, the frame used in the plane-plane case
given some 4-point configuration x (R12).
*/

#ifndef C_PLANES_H
#define C_PLANES_H

#include "SETTINGS.h"
using namespace std;

class C_PLANES
{
    public:
    C_PLANES(){}
    virtual ~C_PLANES() {};

    virtual string name() const
    {
        return string("Collision Frame for vf Plane Distance");
    }

    // utilities for alpha, beta, and gamma derivatives!!
    virtual VECTOR12 dal_dx(const VECTOR12& x) const{
        vector<VECTOR3> v;
        v.resize(4);
        for (int i = 0; i < 4; i++)
        {
            v[i][0] = x[i * 3];
            v[i][1] = x[i * 3 + 1];
            v[i][2] = x[i * 3 + 2];
        }

        VECTOR3 e0, e1, e2;

        e0 = v[1] - v[2];
        e1 = v[3] - v[2];

        VECTOR3 e0h = e0.normalized();
        VECTOR12 da_dx = VECTOR12::Zero();
        da_dx.block<3,1>(3,0) = e1 - 2*e0h*e0h.dot(e1);
        da_dx.block<3,1>(6,0) = 2*e0h*e0h.dot(e1) - e1 - e0;
        da_dx.block<3,1>(9,0) = e0;
        da_dx *= 1/e0.squaredNorm();

        return da_dx;
    }
    virtual VECTOR12 dbet_dx(const VECTOR12& x) const{
        vector<VECTOR3> v;
        v.resize(4);
        for (int i = 0; i < 4; i++)
        {
            v[i][0] = x[i * 3];
            v[i][1] = x[i * 3 + 1];
            v[i][2] = x[i * 3 + 2];
        }

        VECTOR3 e0, e1, e2;

        e0 = v[1] - v[2];
        e1 = v[3] - v[2];
        e2 = v[0] - v[2];

        const REAL e0d0 = e0.dot(e0);
        const REAL e0d1 = e0.dot(e1);
        const REAL e0d2 = e0.dot(e2);
        const REAL e1d1 = e1.dot(e1);
        const REAL e1d2 = e1.dot(e2);
        const REAL e0x12 = (e0.cross(e1)).squaredNorm();

        const REAL bet = (e1d1 * e0d2 - e1d2 * e0d1)/(e0d0 * e1d1 - e0d1 * e0d1);
        
        const VECTOR3 dbdv0 = e1d1*e0 - e0d1*e1;
        const VECTOR3 dbdv1 = 2*bet*(e0d1*e1 - e1d1*e0) - e1d2*e1 + e1d1 *e2;
        const VECTOR3 dbdv2 = 2*bet*(e1d1*e0 + e0d0*e1 - e0d1*(e0+e1)) + e1d2*(e0+e1) - 2*e0d2*e1 - e1d1*(e0+e2) + e0d1*(e1+e2);
        const VECTOR3 dbdv3 = 2*bet*(e0d1*e0 - e0d0*e1) - e1d2*e0 + 2*e0d2*e1 - e0d1*e2;

        VECTOR12 dbdx = VECTOR12::Zero();
        dbdx.block<3,1>(0,0) = dbdv0;
        dbdx.block<3,1>(3,0) = dbdv1;
        dbdx.block<3,1>(6,0) = dbdv2;
        dbdx.block<3,1>(9,0) = dbdv3;

        return dbdx/e0x12;
    }
    virtual VECTOR12 dgam_dx(const VECTOR12& x) const{
        vector<VECTOR3> v;
        v.resize(4);
        for (int i = 0; i < 4; i++)
        {
            v[i][0] = x[i * 3];
            v[i][1] = x[i * 3 + 1];
            v[i][2] = x[i * 3 + 2];
        }

        VECTOR3 e0, e1, e2;

        e0 = v[1] - v[2];
        e1 = v[3] - v[2];
        e2 = v[0] - v[2];

        const REAL e0d0 = e0.dot(e0);
        const REAL e0d1 = e0.dot(e1);
        const REAL e0d2 = e0.dot(e2);
        const REAL e1d1 = e1.dot(e1);
        const REAL e1d2 = e1.dot(e2);
        const REAL e0x12 = (e0.cross(e1)).squaredNorm();

        const REAL gam = (e0d0 * e1d2 - e0d2 * e0d1)/(e0d0 * e1d1 - e0d1 * e0d1);

        const VECTOR3 dgamdv0 = -e0d1*e0 + e0d0*e1;
        const VECTOR3 dgamdv1 = 2*gam*(e0d1*e1 - e1d1*e0) + 2*e1d2*e0 - e0d2*e1 - e0d1*e2;
        const VECTOR3 dgamdv2 = 2*gam*(e1d1*e0 - e0d1*(e0+e1) + e0d0*e1) - 2*e1d2*e0 + e0d2*(e0+e1) + e0d1*(e0+e2) - e0d0*(e1+e2);
        const VECTOR3 dgamdv3 = 2*gam*(e0d1*e0 - e0d0*e1) - e0d2*e0 + e0d0*e2;

        VECTOR12 dgamdx = VECTOR12::Zero();
        dgamdx.block<3,1>(0,0) = dgamdv0;
        dgamdx.block<3,1>(3,0) = dgamdv1;
        dgamdx.block<3,1>(6,0) = dgamdv2;
        dgamdx.block<3,1>(9,0) = dgamdv3;

        return dgamdx/e0x12;
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

        e0 = v[1] - v[2];
        e1 = v[3] - v[2];
        e2 = v[0] - v[2];

        const REAL e0d0 = e0.dot(e0);
        const REAL e0d1 = e0.dot(e1);
        const REAL e0d2 = e0.dot(e2);
        const REAL e1d1 = e1.dot(e1);
        const REAL e1d2 = e1.dot(e2);

        const REAL alph = e0d1/e0d0;
        const REAL bet = (e1d1 * e0d2 - e1d2 * e0d1)/(e0d0 * e1d1 - e0d1 * e0d1);
        const REAL gam = (e0d0 * e1d2 - e0d2 * e0d1)/(e0d0 * e1d1 - e0d1 * e0d1);

        c.block<3,1>(0,0) = e0;
        c.block<3,1>(3,0) = e1 - alph*e0;
        c.block<3,1>(6,0) = e2 - bet*e0 - gam*e1;

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

        e0 = v[1] - v[2];
        e1 = v[3] - v[2];
        e2 = v[0] - v[2];

        const REAL e0d0 = e0.dot(e0);
        const REAL e0d1 = e0.dot(e1);
        const REAL e0d2 = e0.dot(e2);
        const REAL e1d1 = e1.dot(e1);
        const REAL e1d2 = e1.dot(e2);

        const REAL alph = e0d1/e0d0;
        const REAL bet = (e1d1 * e0d2 - e1d2 * e0d1)/(e0d0 * e1d1 - e0d1 * e0d1);
        const REAL gam = (e0d0 * e1d2 - e0d2 * e0d1)/(e0d0 * e1d1 - e0d1 * e0d1);

        MATRIX9x12 partFrame;
        MATRIX3 id = MATRIX3::Identity();
        partFrame.setZero();
        partFrame.block<3,3>(0,3) = id;
        partFrame.block<3,3>(0,6) = -id;
        partFrame.block<3,3>(3,3) = -alph*id;
        partFrame.block<3,3>(3,6) = (alph - 1)*id;
        partFrame.block<3,3>(3,9) = id;
        partFrame.block<3,3>(6,0) = id;
        partFrame.block<3,3>(6,3) = -bet*id;
        partFrame.block<3,3>(6,6) = (bet + gam - 1)*id;
        partFrame.block<3,3>(6,9) = -gam*id;

        return partFrame;
    }

    // A1 term is nonzero with Hu and zero with Hl
    virtual bool has_del(){
        return true;
    }

    virtual MATRIX9x12 del_dCdx_s(const VECTOR12& x) const
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

        e0 = v[1] - v[2];
        e1 = v[3] - v[2];

        VECTOR12 dadx, dbdx, dgdx;
        dadx = dal_dx(x);
        dbdx = dbet_dx(x);
        dgdx = dgam_dx(x);

        MATRIX9x12 partFrame;
        partFrame.setZero();
        partFrame.block<3,3>(3,3) = -e0*dadx.block<3,1>(3,0).transpose();
        partFrame.block<3,3>(3,6) = -e0*dadx.block<3,1>(6,0).transpose();
        partFrame.block<3,3>(3,9) = -e0*dadx.block<3,1>(9,0).transpose();
        partFrame.block<3,3>(6,0) = -e0*dbdx.block<3,1>(0,0).transpose() - e1*dgdx.block<3,1>(0,0).transpose();
        partFrame.block<3,3>(6,3) = -e0*dbdx.block<3,1>(3,0).transpose() - e1*dgdx.block<3,1>(3,0).transpose();
        partFrame.block<3,3>(6,6) = -e0*dbdx.block<3,1>(6,0).transpose() - e1*dgdx.block<3,1>(6,0).transpose();
        partFrame.block<3,3>(6,9) = -e0*dbdx.block<3,1>(9,0).transpose() - e1*dgdx.block<3,1>(9,0).transpose();

        return partFrame;
    }

};

#endif
