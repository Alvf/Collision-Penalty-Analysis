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
#include "TET_WISP_MESH.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Volume/ANISOTROPIC_ARAP.h"
#include "Hyperelastic/Volume/ANISOTROPIC_DIRICHLET.h"
#include "Hyperelastic/Volume/NEGATIVE_DIRICHLET.h"
#include "Hyperelastic/Volume/POWER_DIRICHLET.h"
#include "Hyperelastic/Volume/SNH.h"
#include "Hyperelastic/Volume/TET_STRAND_TWIST.h"
#include "util/MATRIX_UTIL.h"
#include "util/TIMER.h"
#include <iostream>
#include "LINE_INTERSECT.h"
#include <float.h>
#include "Hyperelastic/Volume/ARAP.h"

#define USE_F_TWISTS 1 // Use Alvin's edge twisting formulation? (1 for yes) 

namespace HOBAK {

using namespace std;

///////////////////////////////////////////////////////////////////////
// accepts a vector of individual strands
///////////////////////////////////////////////////////////////////////
TET_WISP_MESH::TET_WISP_MESH(const vector<VECTOR3>& restVertices,
                         const vector<vector<int> >& strandIndices,
                         const REAL& E,        // Young's modulus
                         const REAL& G,       // shear modulus
                         const REAL& density,
                         const REAL& baseRadius,
                         const REAL& tipRadius) :
  TET_STRAND_MESH(restVertices, strandIndices, E, G, density, 0.02, 0.02)
{
  _baseRadius = baseRadius;
  _tipRadius = tipRadius;
  
  //_collisionEps = 0.1;
  //_collisionEps = 2.0 * _tipRadius; // this one given stable results with uniform radii
  //_collisionEps = 0.1 * _tipRadius;

  initializeWisps();

  //cout << " Collision eps:       " << _collisionEps << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TET_WISP_MESH::~TET_WISP_MESH()
{
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
const REAL TET_WISP_MESH::computeTetVolume(const int index) const
{
  vector<VECTOR3> vertices(4);
  const VECTOR4I& tet = _tets[index];
  const VECTOR3& v0 = _vertices[tet[0]];
  const VECTOR3& v1 = _vertices[tet[1]];
  const VECTOR3& v2 = _vertices[tet[2]];
  const VECTOR3& v3 = _vertices[tet[3]];

  const VECTOR3 diff1 = v1 - v0;
  const VECTOR3 diff2 = v2 - v0;
  const VECTOR3 diff3 = v3 - v0;
  return diff3.dot((diff1).cross(diff2)) / 6.0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_WISP_MESH::initializeWisps()
{
  TIMER functionTimer(__FUNCTION__);

  // compute the radii at each vertex
  _vertexRadii.resize(_restVertices.size());
  _vertexRadii.setZero();

  // get the total length of each strand, so we can then interpolate
  _strandLengths.resize(_totalStrands);

  for (unsigned int i = 0; i < _strandIndices.size(); i++)
  {
    const vector<int>& strand = _strandIndices[i];

    // get the total length of the strand
    REAL strandLength = 0;
    for (unsigned int j = 0; j < strand.size() - 1; j++)
    {
      const VECTOR3& v0 = _restVertices[strand[j]];
      const VECTOR3& v1 = _restVertices[strand[j + 1]];

      strandLength += (v0 - v1).norm();
    }
    _strandLengths[i] = strandLength;

    // starting from the base, interpolate from baseRadius to tipRadius
    REAL seenLength = 0;
    for (unsigned int j = 0; j < strand.size(); j++)
    {
      REAL lerp = ((strandLength - seenLength) / strandLength) * _baseRadius + (seenLength / strandLength) * _tipRadius;
      _vertexRadii[strand[j]] = lerp;

      if (j != strand.size() - 1)
      {
        const VECTOR3& v0 = _restVertices[strand[j]];
        const VECTOR3& v1 = _restVertices[strand[j + 1]];
        seenLength += (v0 - v1).norm();
      }
    }
  }

  // propagate vertex radii to edges
  _edgeRadii.resize(_totalEdges);
  _edgeRadii.setZero();
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    // get the edge
    const VECTOR2I edge = _edgeIndices[i];

    // get the radii of each vertex
    const REAL r0 = _vertexRadii[edge[0]];
    const REAL r1 = _vertexRadii[edge[1]];

    _edgeRadii[i] = 0.5 * (r0 + r1);
  }

  // delete previous materials
  _materials.clear();

  // let's make some new ones
  _materials.resize(_totalTets);
  for (int x = 0; x < _totalTets; x++)
  {
    const VECTOR4I tet = _tets[x];

    assert(computeTetVolume[x] > 0.0);

    // TODO: all the vertex lookups
    const REAL r0 = _vertexRadii[tet[0]];
    const REAL r1 = _vertexRadii[tet[1]];
    const REAL r2 = _vertexRadii[tet[2]];
    const REAL r3 = _vertexRadii[tet[3]];

    const REAL r23 = 0.5 * (r2 + r3);
 
    const REAL tetE = _E;
#if 1
    const REAL ks = tetE * M_PI * r23 * r23;
    const REAL kb = tetE * M_PI * (r23 * r23) * 0.25;
    const REAL kt = ((M_PI / 4.0) * _G) * r23 * r23 *
                    (r23 * r23 + r23 * r23);

    /*
    static bool first = true;
    if (first)
    {
      cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
      cout << " ks: " << ks << endl;
      cout << " kb: " << kb << endl;
      cout << " kt: " << kt << endl;
      first = false;
    }
    */
#else
    const REAL ks = 20;
    const REAL kb = 20;
    const REAL kt = 20;
#endif

    // stretching energy for edge near the tip is always added
    _materials[x].stretchingEnergies()[2] = new STRAND::QUADRATIC_STRETCHING(ks);
  
    // add the bending energy 
    const VECTOR3& v0 = _restVertices[tet[0]];
    const VECTOR3& v1 = _restVertices[tet[1]];
    const VECTOR3& v2 = _restVertices[tet[2]];
    const VECTOR3& v3 = _restVertices[tet[3]];
    MATRIX3x2 E;
    E.col(0) = v0 - v1;
    E.col(1) = v2 - v1;
    const REAL theta0 = STRAND::ISOTROPIC_BENDING::angle(E);
    _materials[x].bendingEnergies()[0] = new STRAND::QUADRATIC_UNIT_BENDING(kb, theta0);
    //cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    //cout.precision(16);
    //cout << " bending theta0: " << theta0 << endl;

    // add twist energy for base triangle (Are we using _twistFs?)
    if (x==0) cout << "Using Alvin's edge twist" << endl;
    const VECTOR3& e0 = v0 - v1;
    const VECTOR3& e1 = v2 - v1;
    const VECTOR3& e2 = v3 - v2;
    const VECTOR3& e0Proj = e0 - e0.dot(e1)/e1.squaredNorm() * e1;
    const VECTOR3& e2Proj = e2 - e2.dot(e1)/e1.squaredNorm() * e1;
    const REAL sgn = e0Proj.dot(e1.cross(e2Proj)) < 0 ? -1.0 : 1.0;

    // TK: need to snap this to 1.0 otherwise, it could create a NaN phi0
    REAL projection = e0Proj.normalized().dot(e2Proj.normalized());
    projection = (projection > 1.0) ? 1.0 : projection;
    assert(projection <= 1.0);
    const REAL& phi0 = sgn*acos(projection);
    //cout << " twist phi0: " << phi0 << endl;

    _materials[x].volumeEnergies().push_back(new VOLUME::TET_STRAND_TWIST(kt, phi0));

    // if it's the first tet in a strand, stretching needs to be handled here too
    if (_isStrandBegin[tet[0]])
    {
      const REAL r01 = 0.5 * (r0 + r1);
      const REAL r12 = 0.5 * (r1 + r2);
      const REAL ks01 = _E * M_PI * r01 * r01;
      const REAL ks12 = _E * M_PI * r12 * r12;
 
      _materials[x].stretchingEnergies()[0] = new STRAND::QUADRATIC_STRETCHING(ks01);
      _materials[x].stretchingEnergies()[1] = new STRAND::QUADRATIC_STRETCHING(ks12);
      
      // DEBUG: setting stretching constant monolithiically 
      //_materials[x].stretchingEnergies()[0] = new STRAND::QUADRATIC_STRETCHING(ks);
      //_materials[x].stretchingEnergies()[1] = new STRAND::QUADRATIC_STRETCHING(ks);
    }

    // if it's the last tet in a strand, bending needs to be handled here too
    if (_isStrandEnd[tet[3]])
    {
      const VECTOR3& v3 = _restVertices[tet[3]];
      MATRIX3x2 tipE;
      tipE.col(0) = v1 - v2;
      tipE.col(1) = v3 - v2;
      const REAL tipTheta0 = STRAND::ISOTROPIC_BENDING::angle(tipE);
      _materials[x].bendingEnergies()[1] = new STRAND::QUADRATIC_UNIT_BENDING(kb, tipTheta0);
    }
  }
  _perTetForces.resize(_totalTets);

  computeWispMasses();
  computeStrandIDs();
}
 
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_WISP_MESH::computeStrandIDs()
{
  _edgeStrandID.resize(_totalEdges);
  _edgeStrandID.setZero();

  for (unsigned int i = 0; i < _strandIndices.size(); i++)
  {
    const int strandID = i;
    const vector<int>& strand = _strandIndices[i];

    // find the edge associated with the vertices
    for (unsigned int j = 0; j < strand.size() - 1; j++)
    {
      pair<int,int> edge(strand[j], strand[j+1]);

      // edge index must exist
      auto iter = _edgeHash.find(edge);
      assert(iter != _edgeHash.end());

      // record the strand it's a member of
      _edgeStrandID[iter->second] = strandID;
    }
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_WISP_MESH::computeWispMasses()
{
  assert(_vertexMasses.size() == _totalVertices);

  _vertexMasses.setZero();
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    const REAL& r = _vertexRadii[i];

    _vertexMasses[edge[0]] += 0.5 * _density * M_PI * r * r * _edgeLengths[i];
    _vertexMasses[edge[1]] += 0.5 * _density * M_PI * r * r * _edgeLengths[i];
  }
  
  //for (int x = 0; x < _totalVertices; x++)
  //  cout << " Mass  " << x << ": " << _vertexMasses[x] << endl;
}

}
