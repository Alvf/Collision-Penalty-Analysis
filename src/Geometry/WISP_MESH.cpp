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
#include "WISP_MESH.h"
#include "TIMER.h"
#include <iostream>

using namespace std;

namespace HOBAK {

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
WISP_MESH::WISP_MESH(const vector<VECTOR3>& restVertices,
                     const vector<vector<int> >& strandIndices,
                     const REAL& E,
                     const REAL& G,
                     const REAL& density,
                     const REAL& baseRadius,
                     const REAL& tipRadius) :
  // just send some defaults for density and radius, since they'll get
  // over-written
  STRAND_MESH(restVertices, strandIndices, E, G, density, 0.02, 0.02)
  //STRAND_MESH(restVertices, strandIndices, E, G, 0.0132, 0.02, 0.02)
{
  _baseRadius = baseRadius;
  _tipRadius = tipRadius;
  initializeWisps();
  computeWispMasses();
}

///////////////////////////////////////////////////////////////////////
// generic initialization across multiple constructors
///////////////////////////////////////////////////////////////////////
void WISP_MESH::initializeWisps()
{
  // based on photo of A.M.'s hair, in centimeters
  //const REAL baseRadius = 0.23;
  //const REAL tipRadius = baseRadius / 5.0;
  
  // based on photo of store-bought hair
  //const REAL baseRadius = 1.0 * 0.5;
  //const REAL tipRadius = 0.2 * 0.5;

  // compute the radii at each vertex
  _wispRadii.resize(_restVertices.size());
  _wispRadii.setZero();

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
    //cout << " strand length: " << strandLength << endl;

    // starting from the base, interpolate from baseRadius to tipRadius
    REAL seenLength = 0;
    for (unsigned int j = 0; j < strand.size(); j++)
    {
      REAL lerp = ((strandLength - seenLength) / strandLength) * _baseRadius + (seenLength / strandLength) * _tipRadius;
      _wispRadii[strand[j]] = lerp;

      if (j != strand.size() - 1)
      {
        const VECTOR3& v0 = _restVertices[strand[j]];
        const VECTOR3& v1 = _restVertices[strand[j + 1]];
        seenLength += (v0 - v1).norm();
        //cout << " seen length: " << seenLength << endl;
      }
    }
  }

  //cout << " wisp radii: " << _wispRadii.transpose() << endl;

  // based on the radius, compute the stretching constant
  _kStretch.resize(_totalEdges);
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const VECTOR2I edge = _edgeIndices[x];

    const REAL& r0 = _wispRadii[edge[0]];
    const REAL& r1 = _wispRadii[edge[1]];

    const REAL r = 0.5 * (r0 + r1);
    _kStretch[x] = _E * M_PI * r * r;
    cout << " k stretch " << x << ": " << _kStretch[x] << endl;
  }

  // based on the radius, compute the bending and twisting constant
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const int middleIndex = _bendVertices[x][1];
    const REAL r = _wispRadii[middleIndex];
    const REAL rSq = r * r;

    // bending constant
    MATRIX2 B;
    B.setZero();
    B(0,0) = _E * M_PI * (rSq) * 0.25;
    B(1,1) = _E * M_PI * (rSq) * 0.25;
    _Bs[x] = B;
    cout << " k bend " << x << ": " << B(0,0) << endl;

    // twisting constant
    //const REAL kt = ((1.0 / 8.0) * _E / (1.0 + _nu )) * rSq *
    const REAL kt = ((1.0 / 8.0) * _E / (1.0 + 0.45)) * rSq *
                    (rSq + rSq);
    //const REAL kt = (M_PI * 0.25) * _G * rSq * (rSq + rSq);
    _kts[x] = kt;
  }
  for (unsigned int x = 0; x < _totalBends; x++)
    cout << " k twist " << x << ": " << _kts[x] << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
WISP_MESH::~WISP_MESH()
{
}

///////////////////////////////////////////////////////////////////////
// stretching energy
///////////////////////////////////////////////////////////////////////
REAL WISP_MESH::computeStretchingEnergy(const STRAND::STRETCHING& stretching) const
{
  cout << " Stretching elements: " << endl;
  REAL totalEnergy = 0;
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const VECTOR2I& edge = _edgeIndices[x];
    const REAL k = _kStretch[x];

    // compute deformation gradient
    const VECTOR3 diff = (_vertices[edge[1]] - _vertices[edge[0]]);
    const VECTOR3 f = diff * _dmInvs[x];
    const REAL energy = k * stretching.psi(f);

    cout << x << ": " << energy << " dmInv: " << _dmInvs[x] << " diff.norm() " << diff.norm() << " f.norm(): " << f.norm() << endl;
    totalEnergy += _restEdgeLengths[x] * energy;
  }
  return totalEnergy;
}

///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the stretching force
///////////////////////////////////////////////////////////////////////
VECTOR WISP_MESH::computeStretchingForces(const STRAND::STRETCHING& stretching) const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR6> perEdgeForces(_totalEdges);
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const VECTOR2I edge = _edgeIndices[x];
    const REAL k = _kStretch[x];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    perEdgeForces[x] = -k * _restEdgeLengths[x] * stretching.spatialGradient(positions, 1.0 / _restEdgeLengths[x]);
  }

  return buildPerEdgeVector(perEdgeForces);
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX WISP_MESH::computeStretchingHessian(const STRAND::STRETCHING& stretching) const
{
  vector<MATRIX6> perEdgeHessians(_totalEdges);
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];
    const REAL k = _kStretch[i];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialHessian(positions, 1.0 / _restEdgeLengths[i]);
    perEdgeHessians[i] = -k * H;
  }
  
  return buildPerEdgeMatrix(perEdgeHessians);
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX WISP_MESH::computeStretchingClampedHessian(const STRAND::STRETCHING& stretching) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX6> perEdgeHessians(_totalEdges);
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];
    const REAL k = _kStretch[i];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialClampedHessian(positions, 1.0 / _restEdgeLengths[i]);
    perEdgeHessians[i] = -k * clampEigenvalues(H);
  }

  return buildPerEdgeMatrix(perEdgeHessians);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void WISP_MESH::computeWispMasses()
{
  assert(_vertexMasses.size() == _totalVertices);

  _vertexMasses.setZero();

  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    const REAL& r = _wispRadii[i];

    _vertexMasses[edge[0]] += 0.5 * _density * M_PI * r * r * _edgeLengths[i];
    _vertexMasses[edge[1]] += 0.5 * _density * M_PI * r * r * _edgeLengths[i];
  }

  for (int x = 0; x < _totalVertices; x++)
    cout << " Mass  " << x << ": " << _vertexMasses[x] << endl;
}

}
