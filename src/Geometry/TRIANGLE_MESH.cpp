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
#include "TRIANGLE_MESH.h"
#include "util/MATRIX_UTIL.h"
#include "Collision/COLLISION_UTIL.h"
#include "LINE_INTERSECT.h"
#include "util/TIMER.h"
#include <iostream>
#include <float.h> // need for FLT_MAX

#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"

#define VERY_VERBOSE 0

// Alvin's project: toggle usage of collision energy: 
// 0 for geometry-based sqrt(vf)
// 1 for F-based sqrt(vf)
// 2 for mcAdam's
// 3 for vf
// 4 for F-based vf
// 5 for Adaptive spring constant (F-based wrong for now)
#define FABRIC_USE_COLLISION 0

// Alvin's project: toggling usage of edge collision energies:
// 0 for geometry-based edge-edge
// 1 for F-based edge-edge
// 2 for F-based sqrt(ee)
// 3 for geometry-based sqrtee
// 4 for geometry-bsaed hybrid
#define FABRIC_USE_EDGE 0

namespace HOBAK {

using namespace std;

TRIANGLE_MESH::TRIANGLE_MESH(const vector<VECTOR3>& restVertices, 
                             const vector<VECTOR3I>& triangles) :
    _vertices(restVertices),
    _verticesOld(restVertices),
    _restVertices(restVertices),
    _triangles(triangles)
{
  _DmInvs = computeDmInvs();
  _pFpxs = computePFpxs();

  const int totalTriangles = _triangles.size();
  _Fs.resize(totalTriangles);
  _Us.resize(totalTriangles);
  _Sigmas.resize(totalTriangles);
  _Vs.resize(totalTriangles);

  computeEdges();
  computeSurfaceAreas();
  computeTriangleNeighbors();
  computeEdgeTriangleNeighbors();

  computeFlaps();

  // two centimeters, one seems to get into trouble without CCD
  _collisionEps = 0.02;

  // store which surface vertices are within the one rings of
  // each other
  computeSurfaceVertexOneRings();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TRIANGLE_MESH::~TRIANGLE_MESH()
{
}

///////////////////////////////////////////////////////////////////////
// compute material inverses for deformation gradient
///////////////////////////////////////////////////////////////////////
vector<MATRIX2> TRIANGLE_MESH::computeDmInvs()
{
  vector<MATRIX2> DmInvs(_triangles.size());

  for (unsigned int f = 0; f < _triangles.size(); f++)
  {
    const VECTOR3I& face = _triangles[f];
    const VECTOR3& v0 = _restVertices[face[0]];
    const VECTOR3& v1 = _restVertices[face[1]];
    const VECTOR3& v2 = _restVertices[face[2]];

    MATRIX3x2 Dm;
    Dm.col(0) = v1 - v0;
    Dm.col(1) = v2 - v0;

    // compute the QR factorization
    MATRIX3x2 Qm;
    Qm.col(0) = Dm.col(0).normalized();
    Qm.col(1) = (Dm.col(1) - Dm.col(1).dot(Qm.col(0)) * Qm.col(0)).normalized();

    MATRIX2 Rm;
    Rm = Qm.transpose() * Dm;

    MATRIX2 RmInv;
    const double eps = 1e-8;

    if (fabs(Rm.determinant()) > eps)
      RmInv = Rm.inverse();
    else
    {
      Eigen::JacobiSVD<MATRIX2> svd(Rm, Eigen::ComputeThinU | Eigen::ComputeThinV);
      VECTOR2 singularInvs;
      singularInvs.setZero();
      for (unsigned int x = 0; x < 2; x++)
        if (svd.singularValues()[x] > 1e-8)
          singularInvs[x] = 1.0 / svd.singularValues()[x];

      RmInv = svd.matrixV() * singularInvs.asDiagonal() * svd.matrixU().transpose();
    }
    MATRIX2 DmInv = RmInv;
    DmInvs[f] = DmInv;
  }
  assert(DmInvs.size() == _triangles.size());
  return DmInvs;
}


//////////////////////////////////////////////////////////////////////////////
// derivatives with respect to shape functions
//////////////////////////////////////////////////////////////////////////////
MATRIX3x2 TRIANGLE_MESH::computeDshape(const int i)
{
  MATRIX3x2 dShape;
  dShape.setZero();

  if (i < 3)
  {
    dShape.row(i)[0] = -1;
    dShape.row(i)[1] = -1;
  }
  else if (i < 6)
    dShape.row(i - 3)[0] = 1;
  else
    dShape.row(i - 6)[1] = 1;

  return dShape;
}

///////////////////////////////////////////////////////////////////////
// compute change-of-basis from deformation gradient F to positions, x
///////////////////////////////////////////////////////////////////////
vector<MATRIX6x9> TRIANGLE_MESH::computePFpxs()
{
  vector<MATRIX6x9> pFpxs(_triangles.size());
  pFpxs.clear();
  for (unsigned int f = 0; f < _triangles.size(); f++)
  {
    MATRIX6x9 pFpu;

    MATRIX2 DmInv = _DmInvs[f];
    for (int i = 0; i < 9; i++)
    {
      MATRIX3x2 pDspu = computeDshape(i);
      MATRIX3x2 pFpuColumn = pDspu * DmInv;
      pFpu.col(i) = flatten(pFpuColumn);
    }
    pFpxs.push_back(pFpu);
  }

  return pFpxs;
}

///////////////////////////////////////////////////////////////////////
// for each edge, what're the indices of the neighboring triangles?
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdgeTriangleNeighbors()
{
  // translate the VEC2I into a index
  map<pair<int,int>, int> edgeToIndex;
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    pair<int,int> toHash;
    toHash.first  = _edges[x][0];
    toHash.second = _edges[x][1];

    if (toHash.first > toHash.second)
    {
      int temp = toHash.first;
      toHash.first = toHash.second;
      toHash.second = temp;
    }
    edgeToIndex[toHash] = x;
  }

  // look up the edges of each surface face, tabulate the adjacent
  // triangles
  vector<vector<int> > faceHash(_edges.size());
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I t = _triangles[i];

    // store each edge as pair
    pair<int, int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // store the edge
      assert(edgeToIndex.find(edge) != edgeToIndex.end());
      int edgeIndex = edgeToIndex[edge];
      faceHash[edgeIndex].push_back(i);
    }
  }

  // store the final results
  _edgeTriangleNeighbors.resize(_edges.size());
  for (unsigned int i = 0; i < _edges.size(); i++)
  {
    _edgeTriangleNeighbors[i][0] = -1;
    _edgeTriangleNeighbors[i][1] = -1;

    assert(faceHash[i].size() > 0);
    _edgeTriangleNeighbors[i][0] = faceHash[i][0];

    if (faceHash[i].size() == 2)
      _edgeTriangleNeighbors[i][1] = faceHash[i][1];
  }
}

///////////////////////////////////////////////////////////////////////
// for each triangle, what's the index of the neighboring triangles?
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeTriangleNeighbors()
{
  multimap<pair<int, int>, unsigned int> edgeNeighboringTriangles;

  // hash all the edges from each surface triangle
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I t = _triangles[i];

    // store each edge as pair
    pair<int, int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // hash it
      pair<pair<int, int>, unsigned int> hash(edge, i);
      edgeNeighboringTriangles.insert(hash);
    }
  }
 
  // get the other edge that wasn't the current one 
  _triangleNeighbors.clear();
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I t = _triangles[i];

    // store results here
    VECTOR3I neighbors(-1,-1,-1);

    // reconstruct the edge again
    pair<int,int> edge;
    for (unsigned int j = 0; j < 3; j++)
    {
      edge.first  = t[j];
      edge.second = t[(j + 1) % 3];
      // for cloth, we can have naked edges
      bool orphaned = true; 

      // make sure the ordering is consistent
      if (edge.first > edge.second)
      {
        const int temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      // find the matching triangles
      auto range = edgeNeighboringTriangles.equal_range(edge);
      for (auto it = range.first; it != range.second; it++)
      {
        if (it->second != i){
          neighbors[j] = it->second;
          orphaned = false;
        }
      }
      // if the edgeneighboringtriangles only finds the self triangle
      // input -(j+1).
      if (orphaned) neighbors[j] = -(j+1);
    }

    // store the neighbors
    _triangleNeighbors.push_back(neighbors);
  }
}

///////////////////////////////////////////////////////////////////////
// compute surface areas for collision weights
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSurfaceAreas()
{
  // compute the areas
  _triangleAreas.clear();
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    vector<VECTOR3> vertices(3);
    vertices[0] = _restVertices[_triangles[x][0]];
    vertices[1] = _restVertices[_triangles[x][1]];
    vertices[2] = _restVertices[_triangles[x][2]];

    _triangleAreas.push_back(triangleArea(vertices));
  }

  // cache these out for when triangles deform later
  _restTriangleAreas = _triangleAreas;
  
  // compute the one-ring areas
  assert(_vertices.size() != 0);
  _restOneRingAreas.resize(_vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _restOneRingAreas[x] = 0;
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    assert(x < _triangleAreas.size());
    assert(x < _triangles.size());
    const REAL& area = _triangleAreas[x];
    const VECTOR3I& triangle = _triangles[x];

    for (int y = 0; y < 3; y++)
    {
      const int index = triangle[y];
      assert(index < (int)_restOneRingAreas.size());
      _restOneRingAreas[index] += (1.0 / 3.0) * area;
    }
  }

  // build a mapping from edge index pairs to _edges
  map<pair<int, int>, int> edgeHash;
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    pair<int,int> edge(_edges[x][0], _edges[x][1]);
    edgeHash[edge] = x;
  }

  // compute the edge areas
  assert(_edges.size() != 0);
  _restEdgeAreas.resize(_edges.size());
  _restEdgeAreas.setZero();
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    // build each edge
    for (int y = 0; y < 3; y++)
    {
      pair<int,int> edge(_triangles[x][y],
                         _triangles[x][(y + 1) % 3]);

      // swap them to the order the hash expects
      if (edge.first > edge.second)
      {
        const REAL temp = edge.first;
        edge.first = edge.second;
        edge.second = temp;
      }

      const int edgeIndex = edgeHash[edge];
      assert(edgeIndex >= 0);
      assert(edgeIndex < _restEdgeAreas.size());
      _restEdgeAreas[edgeIndex] += _triangleAreas[x] / 3.0;
    }
  }
}

///////////////////////////////////////////////////////////////////////
// find all edges
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdges()
{
  // hash all the edges, so we don't store any repeats
  map<pair<int, int>, bool> edgeHash;
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++)
    {
      const int v0 = _triangles[x][y];
      const int v1 = _triangles[x][(y + 1) % 3];

      // store them in sorted order
      pair<int, int> edge;
      if (v0 > v1)
      {
        edge.first = v1;
        edge.second = v0;
      }
      else
      {
        edge.first = v0;
        edge.second = v1;
      }

      // hash it out
      edgeHash[edge] = true;
    }
  }

  // store all the unique hashes
  _edges.clear();
  for (auto iter = edgeHash.begin(); iter != edgeHash.end(); iter++)
  {
    const pair<int,int> e = iter->first;
    const VECTOR2I edge(e.first, e.second);
    _edges.push_back(edge);
  }

  cout << " Found " << _edges.size() << " edges on the surface " << endl;
}

///////////////////////////////////////////////////////////////////////
// get all the deformation gradients
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeFs()
{
  TIMER functionTimer(__FUNCTION__);
  assert(_Fs.size() == _triangles.size());

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    _Fs[x] = computeF(x);
    assert(!_Fs[x].hasNaN());
  }

  // SVDs are now stale
  _svdsComputed = false;
}

///////////////////////////////////////////////////////////////////////
// get the SVD of all the deformation gradients
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSVDs()
{
  assert(_Us.size() == _triangles.size());
  assert(_Sigmas.size() == _triangles.size());
  assert(_Vs.size() == _triangles.size());
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _triangles.size(); x++)
    svd(_Fs[x], _Us[x], _Sigmas[x], _Vs[x]);

  _svdsComputed = true;
}

///////////////////////////////////////////////////////////////////////
// get deformation gradient
///////////////////////////////////////////////////////////////////////
MATRIX3x2 TRIANGLE_MESH::computeF(const int index) const
{
  const VECTOR3I& t = _triangles[index];
  MATRIX3x2 Ds;
  Ds.col(0) = _vertices[t[1]] - _vertices[t[0]];
  Ds.col(1) = _vertices[t[2]] - _vertices[t[0]];

  return Ds * _DmInvs[index];
}

///////////////////////////////////////////////////////////////////////
// get stretching energy over entire mesh
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::computeStretchingEnergy(const SHELL::STRETCHING& stretching) const
{
  assert(_triangles.size() == _restTriangleAreas.size());

  VECTOR triangleEnergies(_triangles.size());
  for (int index = 0; index < int(_triangles.size()); index++)
  {
    const MATRIX3x2 F = computeF(index);
    triangleEnergies[index] = _restTriangleAreas[index] * stretching.psi(F);
  }

  return triangleEnergies.sum();
}

///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the stretching force
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeStretchingForces(const SHELL::STRETCHING& stretching) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR9> perElementForces(_triangles.size());
  for (unsigned int index = 0; index < _triangles.size(); index++)
  {
    const MATRIX3x2& F = _Fs[index];
    const MATRIX3x2 PK1 = stretching.PK1(F);
    const VECTOR9 forceDensity = _pFpxs[index].transpose() * flatten(PK1);
    const VECTOR9 force = -_restTriangleAreas[index] * forceDensity;
    perElementForces[index] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int index = 0; index < _triangles.size(); index++)
  {
    const VECTOR3I& triangle = _triangles[index];
    const VECTOR9& triangleForce = perElementForces[index];
    for (int x = 0; x < 3; x++)
    {
      unsigned int index = 3 * triangle[x];
      forces[index]     += triangleForce[3 * x];
      forces[index + 1] += triangleForce[3 * x + 1];
      forces[index + 2] += triangleForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// get bending energy over entire mesh
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::computeBendingEnergy(const SHELL::BENDING& bending) const
{
  assert(_triangles.size() == _restTriangleAreas.size());

  VECTOR flapEnergies(_flaps.size());
  for (int index = 0; index < int(_flaps.size()); index++)
  {
    vector<VECTOR3> flap;
    for (unsigned int j = 0; j < 4; j++)
      flap.push_back(_vertices[_flaps[index][j]]);

    const REAL theta = _restThetas[index];
    flapEnergies[index] = _restFlapAreas[index] * bending.psi(flap, theta);
  }

  return flapEnergies.sum();
}

///////////////////////////////////////////////////////////////////////
// Use the material gradient to compute the stretching force
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeBendingForces(const SHELL::BENDING& bending) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR12> perFlapForces(_flaps.size());
  for (unsigned int index = 0; index < _flaps.size(); index++)
  {
    vector<VECTOR3> flap;
    for (unsigned int j = 0; j < 4; j++)
      flap.push_back(_vertices[_flaps[index][j]]);

    const REAL theta = _restThetas[index];
    const VECTOR12 force = -_restFlapAreas[index] * bending.gradient(flap, theta);
    perFlapForces[index] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    const VECTOR4I& flap      = _flaps[i];
    const VECTOR12& flapForce = perFlapForces[i];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * flap[x];

      if (index + 2 >= DOFs)
      {
        cout << " flap index: " << i << endl;
        cout << " vertex index; " << index << endl;
        cout << " DOFs:  " << DOFs << endl;
      }
      assert(index + 2 < DOFs);
      forces[index]     += flapForce[3 * x];
      forces[index + 1] += flapForce[3 * x + 1];
      forces[index + 2] += flapForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeBendingHessian(const SHELL::BENDING& bending) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perFlapHessians(_flaps.size());
  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    vector<VECTOR3> flap;
    for (unsigned int j = 0; j < 4; j++)
      flap.push_back(_vertices[_flaps[i][j]]);

    const REAL theta = _restThetas[i];
    const MATRIX12 hessian  = -_restFlapAreas[i] * bending.hessian(flap, theta);
    perFlapHessians[i] = hessian;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    const VECTOR4I& flap = _flaps[i];
    const MATRIX12& H = perFlapHessians[i];
    for (int y = 0; y < 4; y++)
    {
      const int yVertex = flap[y];
      for (int x = 0; x < 4; x++)
      {
        const int xVertex = flap[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeBendingClampedHessian(const SHELL::BENDING& bending) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perFlapHessians(_flaps.size());
//#pragma omp parallel
//#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    vector<VECTOR3> flap;
    for (unsigned int j = 0; j < 4; j++)
      flap.push_back(_vertices[_flaps[i][j]]);

    const REAL theta = _restThetas[i];
    const MATRIX12 hessian  = -_restFlapAreas[i] * bending.clampedHessian(flap, theta);
    perFlapHessians[i] = hessian;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _flaps.size(); i++)
  {
    const VECTOR4I& flap = _flaps[i];
    const MATRIX12& H = perFlapHessians[i];
    for (int y = 0; y < 4; y++)
    {
      const int yVertex = flap[y];
      for (int x = 0; x < 4; x++)
      {
        const int xVertex = flap[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeStretchingHessian(const SHELL::STRETCHING& stretching) const
{
  vector<MATRIX9> perElementHessians(_triangles.size());
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const MATRIX3x2& F       = _Fs[i];
    const MATRIX6x9& pFpx = _pFpxs[i];
    const MATRIX6 hessian  = -_restTriangleAreas[i] * stretching.hessian(F);
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I& triangle = _triangles[i];
    const MATRIX9& H = perElementHessians[i];
    for (int y = 0; y < 3; y++)
    {
      const int yVertex = triangle[y];
      for (int x = 0; x < 3; x++)
      {
        const int xVertex = triangle[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeStretchingClampedHessian(const SHELL::STRETCHING& stretching) const
{
  vector<MATRIX9> perElementHessians(_triangles.size());
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const MATRIX3x2& F       = _Fs[i];
    const MATRIX6x9& pFpx = _pFpxs[i];
    const MATRIX6 hessian  = -_restTriangleAreas[i] * stretching.clampedHessian(F);
    perElementHessians[i] = (pFpx.transpose() * hessian) * pFpx;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _triangles.size(); i++)
  {
    const VECTOR3I& triangle = _triangles[i];
    const MATRIX9& H = perElementHessians[i];
    for (int y = 0; y < 3; y++)
    {
      const int yVertex = triangle[y];
      for (int x = 0; x < 3; x++)
      {
        const int xVertex = triangle[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setPositions(const VECTOR& positions)
{
  assert(positions.size() == int(_vertices.size()) * 3);

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] = positions[3 * x];
    _vertices[x][1] = positions[3 * x + 1];
    _vertices[x][2] = positions[3 * x + 2];
  }
}

///////////////////////////////////////////////////////////////////////
// get the current displacement in vector form
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::getDisplacement() const
{
  VECTOR delta(_vertices.size() * 3);
  delta.setZero();

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    const VECTOR3 diff = _vertices[x] - _restVertices[x];
    const int x3 = 3 * x;
    delta[x3]     = diff[0];
    delta[x3 + 1] = diff[1];
    delta[x3 + 2] = diff[2];
  }

  return delta;
}

///////////////////////////////////////////////////////////////////////
// set the vertex displacements to these values exactly
// ALSO UPDDATES _verticesOld
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setDisplacement(const VECTOR& delta)
{
  assert(delta.size() == int(_vertices.size()) * 3); 
  std::swap(_vertices,_verticesOld); //Update the old vertices 
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] = _restVertices[x][0] + delta[3 * x];
    _vertices[x][1] = _restVertices[x][1] + delta[3 * x + 1];
    _vertices[x][2] = _restVertices[x][2] + delta[3 * x + 2];
  }
}

///////////////////////////////////////////////////////////////////////
// add the following deltas to the position
// DOES NOT UPDATE _verticesOld
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::addDisplacement(const VECTOR& delta)
{
  assert(delta.size() == int(_vertices.size()) * 3);

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    _vertices[x][0] += delta[3 * x];
    _vertices[x][1] += delta[3 * x + 1];
    _vertices[x][2] += delta[3 * x + 2];
  }
}

///////////////////////////////////////////////////////////////////////
// read in OBJ-style tet mesh file using Tiny OBJ Loader
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::readObjFile(const string& filename, 
                                vector<VECTOR3>& vertices,
                                vector<VECTOR3I>& triangles, const bool addingOn, const VECTOR3& translate)
{
  // load up the OBJ
  cout << " Reading in *.obj file " << filename.c_str() << endl;

  tinyobj::attrib_t attrib;
  std::vector<tinyobj::shape_t> shapes;
  std::vector<tinyobj::material_t> materials;

  std::string warn;
  std::string err;
  const bool triangulate = true;
  bool success = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str(),
                                  NULL, triangulate);

  if (!warn.empty()) {
    cout << "WARNING: " << warn << endl;
  }
  if (!err.empty()) {
    cout << "ERROR: " << err << endl;
  }
  if (!success) {
    cout << " Failed to read " << filename.c_str() << endl;
    return false;
  }

  // erase whatever was in the vectors before
  if(!addingOn){
    vertices.clear();
    triangles.clear();
  }

  int base_offset = vertices.size();

  // store up the vertices
  const vector<float>& vs = attrib.vertices;
  for (unsigned int x = 0; x < (vs.size() / 3); x++)
  {
    const int x3 = 3 * x;
    const VECTOR3 v(vs[x3 + 0] + translate[0], vs[x3 + 1] + translate[1], vs[x3 + 2] + translate[2]);
    vertices.push_back(v);
  }

  // there's at least one shape, right?
  assert(shapes.size() > 0);

  // just pick off the first shape
  size_t index_offset = 0;
  for (size_t f = 0; f < shapes[0].mesh.num_face_vertices.size(); f++)
  {
    VECTOR3I triangle;

    // it's a triangle, right?
    const size_t fnum = shapes[0].mesh.num_face_vertices[f];
    assert(fnum == 3);
    for (size_t v = 0; v < 3; v++) 
    {
      tinyobj::index_t idx = shapes[0].mesh.indices[index_offset + v];
      int index = idx.vertex_index;

      triangle[v] = index + base_offset;
    }
    index_offset += 3;

    triangles.push_back(triangle);
  }
  cout << "total " << vertices.size() << " vertices " << endl;
  cout << " total " << triangles.size() << " faces " << endl;

  return true;
}

///////////////////////////////////////////////////////////////////////
// write out OBJ-style tet mesh file
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::writeObjFile(const string& filename, const TRIANGLE_MESH& triMesh, 
                                 const bool restVertices)
{
  FILE* file = fopen(filename.c_str(), "w");
  cout << " Writing out mesh file: " << filename.c_str() << endl;
 
  if (file == NULL)
  {
    cout << " Failed to open file!" << endl;
    return false;
  }

  const vector<VECTOR3>& vertices = (restVertices) ? triMesh.restVertices() : triMesh.vertices(); 
  const vector<VECTOR3I>& tris = triMesh.triangles(); 

  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    const VECTOR3& v = vertices[x];
    fprintf(file, "v %.17g %.17g %.17g\n", v[0], v[1], v[2]);
  }
  for (unsigned int x = 0; x < tris.size(); x++)
  {
    const VECTOR3I& tri = tris[x];
    fprintf(file, "f %i %i %i\n", tri[0] + 1, tri[1] + 1, tri[2] + 1);
  }

  fclose(file);
  return true;
}

///////////////////////////////////////////////////////////////////////
// normalize vertices so that they're in a unit box, 
// centered at (0.5, 0.5, 0.5)
///////////////////////////////////////////////////////////////////////
vector<VECTOR3> TRIANGLE_MESH::normalizeVertices(const vector<VECTOR3>& vertices)
{
  assert(vertices.size() > 0);
  VECTOR3 mins = vertices[0];
  VECTOR3 maxs = vertices[0];
  for (unsigned int x = 1; x < vertices.size(); x++)
    for (int y = 0; y < 3; y++)
    {
      mins[y] = (mins[y] < vertices[x][y]) ? mins[y] : vertices[x][y];
      maxs[y] = (maxs[y] > vertices[x][y]) ? maxs[y] : vertices[x][y];
    }

  const VECTOR3 lengths = maxs - mins;
  const REAL maxLengthInv = 1.0 / lengths.maxCoeff();

  vector<VECTOR3> normalized = vertices;
  for (unsigned int x = 0; x < vertices.size(); x++)
  {
    normalized[x] -= mins;
    normalized[x] *= maxLengthInv;
    
    normalized[x] += VECTOR3(0.5, 0.5, 0.5);
  }

  return normalized;
}

///////////////////////////////////////////////////////////////////////
// see if the projection of v onto the plane of v0,v1,v2 is inside 
// the triangle formed by v0,v1,v2
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::pointProjectsInsideTriangle(const VECTOR3& v0, const VECTOR3& v1, 
                                           const VECTOR3& v2, const VECTOR3& v)
{
  // get the barycentric coordinates
  const VECTOR3 e1 = v1 - v0;
  const VECTOR3 e2 = v2 - v0;
  const VECTOR3 n = e1.cross(e2);
  const VECTOR3 na = (v2 - v1).cross(v - v1);
  const VECTOR3 nb = (v0 - v2).cross(v - v2);
  const VECTOR3 nc = (v1 - v0).cross(v - v0);
  const VECTOR3 barycentric(n.dot(na) / n.squaredNorm(),
                            n.dot(nb) / n.squaredNorm(),
                            n.dot(nc) / n.squaredNorm());
                            
  const REAL barySum = fabs(barycentric[0]) + fabs(barycentric[1]) + fabs(barycentric[2]);

  // if the point projects to inside the triangle, it should sum to 1
  if (barySum - 1.0 < 1e-8)
    return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute distance between a point and triangle
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::pointTriangleDistance(const VECTOR3& v0, const VECTOR3& v1, 
                                     const VECTOR3& v2, const VECTOR3& v)
{
  // get the barycentric coordinates
  const VECTOR3 e1 = v1 - v0;
  const VECTOR3 e2 = v2 - v0;
  const VECTOR3 n = e1.cross(e2);
  const VECTOR3 na = (v2 - v1).cross(v - v1);
  const VECTOR3 nb = (v0 - v2).cross(v - v2);
  const VECTOR3 nc = (v1 - v0).cross(v - v0);
  const VECTOR3 barycentric(n.dot(na) / n.squaredNorm(),
                            n.dot(nb) / n.squaredNorm(),
                            n.dot(nc) / n.squaredNorm());
                            
  const REAL barySum = fabs(barycentric[0]) + fabs(barycentric[1]) + fabs(barycentric[2]);

  // if the point projects to inside the triangle, it should sum to 1
  if (barySum - 1.0 < 1e-8)
  {
    const VECTOR3 nHat = n / n.norm();
    const REAL normalDistance = (nHat.dot(v - v0));
    return fabs(normalDistance);
  }

  // project onto each edge, find the distance to each edge
  const VECTOR3 e3 = v2 - v1;
  const VECTOR3 ev = v - v0;
  const VECTOR3 ev3 = v - v1;
  const VECTOR3 e1Hat = e1 / e1.norm();
  const VECTOR3 e2Hat = e2 / e2.norm();
  const VECTOR3 e3Hat = e3 / e3.norm();
  VECTOR3 edgeDistances(FLT_MAX, FLT_MAX, FLT_MAX);

  // see if it projects onto the interval of the edge
  // if it doesn't, then the vertex distance will be smaller,
  // so we can skip computing anything
  const REAL e1dot = e1Hat.dot(ev);
  if (e1dot > 0.0 && e1dot < e1.norm())
  {
    const VECTOR3 projected = v0 + e1Hat * e1dot;
    edgeDistances[0] = (v - projected).norm();
  }
  const REAL e2dot = e2Hat.dot(ev);
  if (e2dot > 0.0 && e2dot < e2.norm())
  {
    const VECTOR3 projected = v0 + e2Hat * e2dot;
    edgeDistances[1] = (v - projected).norm();
  }
  const REAL e3dot = e3Hat.dot(ev3);
  if (e3dot > 0.0 && e3dot < e3.norm())
  {
    const VECTOR3 projected = v1 + e3Hat * e3dot;
    edgeDistances[2] = (v - projected).norm();
  }

  // get the distance to each vertex
  const VECTOR3 vertexDistances((v - v0).norm(), 
                                (v - v1).norm(), 
                                (v - v2).norm());

  // get the smallest of both the edge and vertex distances
  const REAL vertexMin = vertexDistances.minCoeff();
  const REAL edgeMin = edgeDistances.minCoeff();

  // return the smallest of those
  return (vertexMin < edgeMin) ? vertexMin : edgeMin;
}

///////////////////////////////////////////////////////////////////////
// see if the vertex is inside the collision cell described in 
// Chapter 11: Collision Processing, [Kim and Eberle 2020]
//
// this seems to work okay, but need stronger debugging checks
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::insideCollisionCell(const int triangleID, const VECTOR3& vertex)
{
  const VECTOR3I& t = _triangles[triangleID];
  vector<VECTOR3> v;
  v.push_back(_vertices[t[0]]);
  v.push_back(_vertices[t[1]]);
  v.push_back(_vertices[t[2]]);
  VECTOR3 n = planeNormal(v);

  // get the normals of the three adjacent faces
  const VECTOR3I& neighbors = _triangleNeighbors[triangleID];
  for (int x = 0; x < 3; x++)
  {
    VECTOR3 nebHat;
    // Just skip if it's boundary edge; distance checks in higher function
    // will make sure nothing crazy happens
    if (neighbors[x] < 0){
      // int j = -neighbors[x] -1;
      // VECTOR3 p0 = (v[(j+1)%3] + v[j])/2;
      // VECTOR3 p1 = v[(j+2)%3];
      // VECTOR3 edge = (v[j+1%3] - v[j]).normalized();
      // nebHat = (p1-p0 -(p1-p0).dot(edge)*edge).normalized();
      continue;
    }
    else{
      const VECTOR3I& tNeighbor = _triangles[neighbors[x]];
      vector<VECTOR3> vNeighbor;
      vNeighbor.push_back(_vertices[tNeighbor[0]]);
      vNeighbor.push_back(_vertices[tNeighbor[1]]);
      vNeighbor.push_back(_vertices[tNeighbor[2]]);
      const VECTOR3 nNeighbor = planeNormal(vNeighbor);
      // do the inside check
      const VECTOR3 ne = (nNeighbor + n).normalized();
      const VECTOR3 eij = v[(x+1) % 3] - v[x];
      const VECTOR3 neb = ne.cross(eij);
      nebHat = neb.normalized();
    }

    const REAL deplane = nebHat.dot(vertex - v[x]);
    if (deplane < 0.0)
      return false;
  }

  // TODO: some face plane compares, and returning if the vertex is above
  // or below the face
  return true;
}

// ///////////////////////////////////////////////////////////////////////
// // compute distance to collision cell wall, where positive means inside
// // and negative means outside
// // does this get used anywhere?
// ///////////////////////////////////////////////////////////////////////
// REAL TRIANGLE_MESH::distanceToCollisionCellWall(const int triangleID, const VECTOR3& vertex)
// {
//   const VECTOR3I& t = _triangles[triangleID];
//   vector<VECTOR3> v;
//   v.push_back(_vertices[t[0]]);
//   v.push_back(_vertices[t[1]]);
//   v.push_back(_vertices[t[2]]);
//   VECTOR3 n = planeNormal(v);

//   REAL smallestDistance = FLT_MAX;
//   // get the normals of the three adjacent faces
//   vector<VECTOR3> nNeighbors;
//   const VECTOR3I& neighbors = _triangleNeighbors[triangleID];
//   for (int x = 0; x < 3; x++)
//   {
//     VECTOR3 nebHat;
//     // inside check is pretty different if considering orphaned or unorphaned edge
//     if (neighbors[x] < 0){
//       int j = -neighbors[x] -1;
//       VECTOR3 p0 = (v[(j+1)%3] + v[j])/2;
//       VECTOR3 p1 = v[(j+2)%3];
//       nebHat = (p1-p0).normalized();
//     }
//     else{
//       const VECTOR3I& tNeighbor = _triangles[neighbors[x]];
//       vector<VECTOR3> vNeighbor;
//       vNeighbor.push_back(_vertices[tNeighbor[0]]);
//       vNeighbor.push_back(_vertices[tNeighbor[1]]);
//       vNeighbor.push_back(_vertices[tNeighbor[2]]);
//       const VECTOR3 nNeighbor = planeNormal(vNeighbor);
//       // do the inside check
//       const VECTOR3 ne = (nNeighbor + n).normalized();
//       const VECTOR3 eij = v[(x+1) % 3] - v[x];
//       const VECTOR3 neb = ne.cross(eij);
//       nebHat = neb.normalized();
//     }

//     const REAL deplane = nebHat.dot(vertex - v[x]);
//     if (fabs(deplane) < smallestDistance)
//     {
//       smallestDistance = fabs(deplane);
//     }
//   }

//   return smallestDistance;
// }


///////////////////////////////////////////////////////////////////////
// find all the vertex-face collision pairs, using the 
// InFaceRegion test from "Collision Processing" chapter of
// "Dynamic Deformables"
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeVertexFaceCollisions()
{
  TIMER functionTimer(__FUNCTION__);

  // if a vertex is part of an inverted triangle, don't have it participate 
  // in a self-collision. That triangle needs to get its house in order 
  // before it starts bossing around a surface face. Not checking for 
  // this causes faces to get horribly tangled in inverted configurations.
  computeInvertedVertices();

  _vertexFaceCollisions.clear();
  const REAL collisionEps = _collisionEps;

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    // if the vertex is involved in an inverted triangle, give up
    if (_invertedVertices[x]) 
      continue;

    const VECTOR3& surfaceVertex = _vertices[x];

    // find the close triangles
    for (unsigned int y = 0; y < _triangles.size(); y++)
    {
      // if the surface triangle is so small the normal could be degenerate, skip it
      if (surfaceTriangleIsDegenerate(y))
        continue;

      const VECTOR3I& t = _triangles[y];

      // if it's an inverted face, move on
      if (_invertedVertices[t[0]] && _invertedVertices[t[1]] && _invertedVertices[t[2]])
        continue;

      // if this triangle is in the one-ring of the current vertex, skip it
      if (t[0] == x || t[1] == x || t[2] == x) continue;
      
      const REAL distance = pointTriangleDistance(_vertices[t[0]], _vertices[t[1]],
                                                  _vertices[t[2]], surfaceVertex);

      if (distance < collisionEps)
      {
        // if the point, projected onto the face's plane, is inside the face,
        // then record the collision now
        if (pointProjectsInsideTriangle(_vertices[t[0]], _vertices[t[1]], 
                                        _vertices[t[2]], surfaceVertex))
        {
          pair<int,int> collision(x, y);
          _vertexFaceCollisions.push_back(collision);
          continue;
        }
        if (insideCollisionCell(y, surfaceVertex))
        {
          pair<int,int> collision(x, y);
          _vertexFaceCollisions.push_back(collision);
        }
      }
    }
  }

if (_vertexFaceCollisions.size() > 0){
  cout << " Found " << _vertexFaceCollisions.size() << " vertex-face collisions " << endl;
}
#if VERY_VERBOSE
  if (_vertexFaceCollisions.size() > 0){
    cout << " Found " << _vertexFaceCollisions.size() << " vertex-face collisions " << endl;
    cout << " pairs: " << endl;
  }
  for (unsigned int x = 0; x < _vertexFaceCollisions.size(); x++)
  {
    const pair<int,int> collision = _vertexFaceCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")" << endl;
  }
#endif
}

///////////////////////////////////////////////////////////////////////
// find vertex-face CCD collisions 
// through comparing _verticesOld with _vertices
///////////////////////////////////////////////////////////////////////
// void TRIANGLE_MESH::computeVertexFaceCCD()
// {
// #ifndef DISABLE_CCD
//   TIMER functionTimer(__FUNCTION__);

//   // if a vertex is part of an inverted triangle, don't have it participate 
//   // in a self-collision. That triangle needs to get its house in order 
//   // before it starts bossing around a surface face. Not checking for 
//   // this causes faces to get horribly tangled in inverted configurations.
//   computeInvertedVertices();

//   _vertexFaceCCDs.clear(); // data structure for storing CCD vf collisions
//   const REAL collisionEps = _collisionEps;

//   for (unsigned int x = 0; x < _vertices.size(); x++)
//   {
//     // if the vertex is involved in an inverted triangle, give up
//     if (_invertedVertices[x]) 
//       continue;

//     // find the close triangles
//     for (unsigned int y = 0; y < _triangles.size(); y++)
//     {
//       // if the surface triangle is so small the normal could be degenerate, skip it
//       if (surfaceTriangleIsDegenerate(y))
//         continue;

//       const VECTOR3I& t = _triangles[y];

//       // if it's an inverted face, move on
//       if (_invertedVertices[t[0]] && _invertedVertices[t[1]] && _invertedVertices[t[2]])
//         continue;

//       // if this triangle is in the one-ring of the current vertex, skip it
//       if (t[0] == x || t[1] == x || t[2] == x) continue;
      
//       // Do the CCD test: If CCD succeeds, save the vertex index and the triangle index
      
//       // The colliding vertex
//       REAL* x01 = _vertices[x].data();
//       REAL* x00 = _verticesOld[x].data();
//       // The triangle
//       REAL* x10 = _verticesOld[t[0]].data();
//       REAL* x11 = _vertices[t[0]].data();
//       REAL* x20 = _verticesOld[t[1]].data();
//       REAL* x21 = _vertices[t[1]].data();
//       REAL* x30 = _verticesOld[t[2]].data();
//       REAL* x31 = _vertices[t[2]].data();

//       REAL tHit;

//       if(_ccd->Vertex_Triangle_CCD(x00,x01,x10,x11,x20,x21,x30,x31,tHit)){
//         pair<REAL,pair<int,int>> collisionPair;
//         collisionPair.second.first = x;
//         collisionPair.second.second = y;
//         collisionPair.first = tHit;
//         _vertexFaceCCDs.push_back(collisionPair);
//       }

//     }
//   }

//   // sort the vf CCDs by time
//   std::sort(_vertexFaceCCDs.begin(), _vertexFaceCCDs.end());

// #if VERY_VERBOSE
//   if (_vertexFaceCCDs.size() > 0){
//     cout << " Found " << _vertexFaceCCDs.size() << " vertex-face CCDs " << endl;
//     cout << " pairs: " << endl;
//   }
//   for (unsigned int x = 0; x < _vertexFaceCCDs.size(); x++)
//   {
//     const pair<REAL, pair<int,int>> collision = _vertexFaceCCDs[x];
//     cout << "(" << collision.first << ", " << collision.second.first << ", " << collision.second.second << ")" << endl;
//   }
// #endif
//   #endif
//   return;
// }

///////////////////////////////////////////////////////////////////////
// find all the edge-edge self collision pairs, using the 
// brute-force tests
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeEdgeEdgeCollisions()
{
  TIMER functionTimer(__FUNCTION__);
  _edgeEdgeCollisions.clear();
  _edgeEdgeIntersections.clear();
  //_edgeEdgeCollisionEps.clear();
  _edgeEdgeCoordinates.clear();
  _edgeEdgeCollisionAreas.clear();

  // build a mapping from edge index pairs to _edges
  map<pair<int, int>, int> edgeHash;
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    pair<int,int> edge(_edges[x][0], _edges[x][1]);
    edgeHash[edge] = x;
  }

  // get the nearest edge to each edge, not including itself
  // and ones where it shares a vertex
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    int closestEdge = -1;
    REAL closestDistance = FLT_MAX;
    VECTOR2 aClosest(-1,-1);
    VECTOR2 bClosest(-1,-1);
    const VECTOR2I outerEdge = _edges[x];
    const VECTOR3& v0 = _vertices[outerEdge[0]];
    const VECTOR3& v1 = _vertices[outerEdge[1]];

    //const unsigned int outerFlat = outerEdge[0] + outerEdge[1] * _edges.size();
    
    // find the closest other edge
    for (unsigned int y = x + 1; y < _edges.size(); y++)
    {
      const VECTOR2I innerEdge = _edges[y];
      // if they share a vertex, skip it
      if ((outerEdge[0] == innerEdge[0]) || (outerEdge[0] == innerEdge[1]) ||
          (outerEdge[1] == innerEdge[0]) || (outerEdge[1] == innerEdge[1]))
        continue;

      const VECTOR3& v2 = _vertices[innerEdge[0]];
      const VECTOR3& v3 = _vertices[innerEdge[1]];

      // call the geometric test
      VECTOR3 innerPoint, outerPoint;
      IntersectLineSegments(v0, v1, v2, v3,
                            outerPoint, innerPoint);  

      const REAL distance = (innerPoint - outerPoint).norm();

      //if (distance < separationDistance)
      //  separationDistance = distance;

      if (distance > closestDistance) continue;

      // get the line interpolation coordinates
      VECTOR2 a,b;
      const VECTOR3 e0 = v1 - v0;
      const VECTOR3 e1 = v3 - v2;

      // this is a little dicey in general, but if the intersection test isn't
      // total garbage, it should still be robust
      a[1] = (outerPoint - v0).norm() / e0.norm();
      a[0] = 1.0 - a[1];
      b[1] = (innerPoint - v2).norm() / e1.norm();
      b[0] = 1.0 - b[1];

      // if it's really close to an end vertex, skip it
      const REAL skipEps = 1e-4;
      if ((a[0] < skipEps) || (a[0] > 1.0 - skipEps)) continue;
      if ((a[1] < skipEps) || (a[1] > 1.0 - skipEps)) continue;
      if ((b[0] < skipEps) || (b[0] > 1.0 - skipEps)) continue;
      if ((b[1] < skipEps) || (b[1] > 1.0 - skipEps)) continue;

      // it's mid-segment, and closest, so remember it
      closestDistance = distance;
      closestEdge = y;

      aClosest = a;
      bClosest = b;
    }

    // if nothing was close, move on
    if (closestEdge == -1) continue;

    // are they within each other's one rings?
    const VECTOR2I innerEdge = _edges[closestEdge];
    bool insideOneRing = false;

    for (int j = 0; j < 2; j++)
    {
      pair<int, int> lookup;
      lookup.first = outerEdge[j];
      for (int i = 0; i < 2; i++)
      {
        lookup.second = innerEdge[i];
        if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end())
          insideOneRing = true;
      }
    }
    if (insideOneRing) continue;

    // if it's within the positive threshold, it's in collision
    if (closestDistance < _collisionEps)
    {
      pair<int,int> collision(x, closestEdge);
      _edgeEdgeCollisions.push_back(collision);
      
      // this was actually set, right?
      assert(aClosest[0] > 0.0 && aClosest[1] > 0.0);
      assert(bClosest[0] > 0.0 && bClosest[1] > 0.0);

      pair<VECTOR2,VECTOR2> coordinate(aClosest, bClosest);
      _edgeEdgeCoordinates.push_back(coordinate);

      // get the areas too
      const VECTOR2I innerEdge = _edges[closestEdge];
      const pair<int,int> outerPair(outerEdge[0], outerEdge[1]);
      const pair<int,int> innerPair(innerEdge[0], innerEdge[1]);
      const REAL xArea = _restEdgeAreas[edgeHash[outerPair]];
      const REAL closestArea = _restEdgeAreas[edgeHash[innerPair]];
      _edgeEdgeCollisionAreas.push_back(xArea + closestArea);

      // find out if they are penetrating
      vector<VECTOR3> edge(2);
      edge[0] = v0;
      edge[1] = v1;

      // get the adjacent triangles of the *other* edge
      VECTOR2I adjacentTriangles = _edgeTriangleNeighbors[edgeHash[innerPair]];

      // build triangle 0
      const VECTOR3I surfaceTriangle0 = _triangles[adjacentTriangles[0]];
      vector<VECTOR3> triangle0;
      triangle0.push_back(_vertices[surfaceTriangle0[0]]);
      triangle0.push_back(_vertices[surfaceTriangle0[1]]);
      triangle0.push_back(_vertices[surfaceTriangle0[2]]);

      // build triangle 1
      vector<VECTOR3> triangle1;
      // if there's another triangle on the other side (this is in case we're looking at cloth)
      // then store that one too
      if (adjacentTriangles[1] != -1)
      {
        const VECTOR3I surfaceTriangle1 = _triangles[adjacentTriangles[1]];
        triangle1.push_back(_vertices[surfaceTriangle1[0]]);
        triangle1.push_back(_vertices[surfaceTriangle1[1]]);
        triangle1.push_back(_vertices[surfaceTriangle1[2]]);
      }

      // see if the edges are already penetrating the opposing faces
      bool penetrating = false;
      if (triangle0.size() > 0) penetrating = faceEdgeIntersection(triangle0, edge);
      if (triangle1.size() > 0) penetrating = penetrating || faceEdgeIntersection(triangle1, edge);

      _edgeEdgeIntersections.push_back(penetrating);

      // TODO: for completeness, should probably test the other edges against the other
      // pair, just in case we're looking at a degenerate case. In general, seems redundant.
    }
  }
  assert(_edgeEdgeCollisions.size() == _edgeEdgeCoordinates.size());

  if (_edgeEdgeCollisions.size() > 0){
    cout << " Found " << _edgeEdgeCollisions.size() << " edge-edge collisions " << endl;
  }
#if VERY_VERBOSE
  if (_edgeEdgeCollisions.size() > 0){
    cout << " Found " << _edgeEdgeCollisions.size() << " edge-edge collisions " << endl;
    cout << " Collision area array size: " << _edgeEdgeCollisionAreas.size() << endl;
    cout << " Intersections array size: " << _edgeEdgeIntersections.size() << endl;
    cout << " pairs: " << endl;
  }
  for (unsigned int x = 0; x < _edgeEdgeCollisions.size(); x++)
  {
    const pair<int,int> collision = _edgeEdgeCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")" << endl;
  }
#endif
}

///////////////////////////////////////////////////////////////////////
// get the triangle area
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::triangleArea(const vector<VECTOR3>& triangle)
{
  const VECTOR3 edge1 = triangle[1] - triangle[0];
  const VECTOR3 edge2 = triangle[2] - triangle[0];
  return 0.5 * edge1.cross(edge2).norm();
}

///////////////////////////////////////////////////////////////////////
// get the normal to a plane, specified by three points
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::planeNormal(const vector<VECTOR3>& plane)
{
  const VECTOR3 edge1 = plane[1] - plane[0];
  const VECTOR3 edge2 = plane[2] - plane[0];
  return edge1.cross(edge2).normalized();
}

///////////////////////////////////////////////////////////////////////
// project point onto plane, specific by three points
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::pointPlaneProjection(const vector<VECTOR3>& plane, const VECTOR3& point)
{
  const VECTOR3 normal = planeNormal(plane);
  return point - (normal.dot(point - plane[0])) * normal;
}


///////////////////////////////////////////////////////////////////////
// based on vertex-face collision pairs, build "collision tets"
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::buildVertexFaceCollisionTets(const VECTOR& velocity)
{
  // clear out old ones (not clear if this is smart or dumb)
  _vertexFaceCollisionTets.clear();
  _vertexFaceCollisionAreas.clear();

  // make a tet for each vertex-face pair
  for (unsigned int x = 0; x < _vertexFaceCollisions.size(); x++)
  {
    const int vertexID = _vertexFaceCollisions[x].first;
    const int faceID = _vertexFaceCollisions[x].second;
    const VECTOR3I& face = _triangles[faceID];
    
    vector<VECTOR3> vs(4);
    for(int i = 0; i < 3; i++){
      vs[i] = _vertices[face[i]];
    }

    // build a tet with the correct vertex ordering
    VECTOR4I tet;
    tet[0] = vertexID;

    // cloth collisions: if face normal is aligned with collision
    // edge, invert face order. Otherwise feed face order normally

    VECTOR3 faceNorm = planeNormal(vs);
    VECTOR3 wTests = _vertices[vertexID] - _vertices[face[1]];
    if(wTests.dot(faceNorm) > 0)
    {
      tet[1] = face[2];
      tet[2] = face[1];
      tet[3] = face[0];
    }
    else
    {
      tet[1] = face[0];
      tet[2] = face[1];
      tet[3] = face[2];
    }
    // get the rest area of the triangle
    vector<VECTOR3> restFace(3);
    restFace[0] = _restVertices[face[0]]; 
    restFace[1] = _restVertices[face[1]];
    restFace[2] = _restVertices[face[2]];
    const REAL restFaceArea = triangleArea(restFace);

    const REAL restVertexArea = _restOneRingAreas[vertexID];

    // store everything
    _vertexFaceCollisionTets.push_back(tet);
    assert(restFaceArea >= 0.0);
    assert(restVertexArea >= 0.0);
    _vertexFaceCollisionAreas.push_back(restFaceArea + restVertexArea);
  }
}

///////////////////////////////////////////////////////////////////////
// compute collision forces using collision tets
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeVertexFaceCollisionForces(const ENERGY_12D& energy) const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perElementForces(_vertexFaceCollisionTets.size());
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    vector<VECTOR3> vs(4);
    for (unsigned int j = 0; j < 4; j++){
      vs[j] = _vertices[_vertexFaceCollisionTets[i][j]];
    }
    const VECTOR12 force = -_vertexFaceCollisionAreas[i] * energy.gradient(vs);
    perElementForces[i] = force;

#if ENABLE_DEBUG_TRAPS
    if (force.hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " NaN in collision tet: " << i << endl;
      for (int j = 0; j < 4; j++)
        cout << " v" << j << ": " << vs[j].transpose() << endl;
      cout << " gradient: " << endl << _vertexFaceEnergy->gradient(vs) << endl;
    }
#endif
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    const VECTOR4I& tet = _vertexFaceCollisionTets[i];
    const VECTOR12& tetForce = perElementForces[i];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * tet[x];
      forces[index]     += tetForce[3 * x];
      forces[index + 1] += tetForce[3 * x + 1];
      forces[index + 2] += tetForce[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision energy using x-based formulation
///////////////////////////////////////////////////////////////////////
REAL TRIANGLE_MESH::computeEdgeEdgeCollisionEnergy(const ENERGY_12D& eeEnergy) const
{
  TIMER functionTimer(__FUNCTION__);

  REAL finalEnergy = 0.0;
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];
    VECTOR3 e0n = (vs[1] - vs[0]).normalized();
    VECTOR3 e1n = (vs[3] - vs[2]).normalized();
    if (abs(e0n.dot(e1n)) > 1 - 1e-2) {//near-parallel guard
      continue;  
    }
    vector<int> extraInfo(1);
    extraInfo[0] = 1;
    REAL len = eeEnergy.getLen(eeEnergy.makeState(vs));
    bool isNeg = (len < 0);
    bool intersect = _edgeEdgeIntersections[i];
    if(isNeg != intersect){
      extraInfo[0] = -1;
    }

    const REAL psi = eeEnergy.psi(vs,extraInfo);
    finalEnergy += _edgeEdgeCollisionAreas[i] * psi;
  }

  return finalEnergy;
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision forces using x-based formulation
///////////////////////////////////////////////////////////////////////
VECTOR TRIANGLE_MESH::computeEdgeEdgeCollisionForces(const ENERGY_12D& edgeEdgeEnergy) const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perElementForces(_edgeEdgeCollisions.size());
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];
    VECTOR3 e0n = (vs[1] - vs[0]).normalized();
    VECTOR3 e1n = (vs[3] - vs[2]).normalized();
    if (abs(e0n.dot(e1n)) > 1 - 1e-2) {//near-parallel guard
      perElementForces[i] = VECTOR12::Zero();
      continue;  
    }

    REAL len = edgeEdgeEnergy.getLen(edgeEdgeEnergy.makeState(vs));
    bool isNeg = (len < 0);
    bool intersect = _edgeEdgeIntersections[i];
    vector<int> extraInfo(1);
    extraInfo[0] = 1;
    if(isNeg != intersect){
      extraInfo[0] = -1;
    }
    REAL area = _edgeEdgeCollisionAreas[i];
    VECTOR12 g = edgeEdgeEnergy.gradient(vs,extraInfo);
    const VECTOR12 force = -area * g;
    perElementForces[i] = force;
  }

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  { 
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];
    VECTOR12 force = perElementForces[i];
    vector<int> vertexIndices(4);
    vertexIndices[0] = edge0[0];
    vertexIndices[1] = edge0[1];
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];
    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * vertexIndices[x];
      assert((int)index < DOFs);
      forces[index]     += force[3 * x];
      forces[index + 1] += force[3 * x + 1];
      forces[index + 2] += force[3 * x + 2];
    }
  }
  
  return forces;
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision Hessians using x-based formulation
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeEdgeEdgeCollisionClampedHessian(const ENERGY_12D& eeEnergy) const
{
  TIMER functionTimer(__FUNCTION__);

  vector<MATRIX12> perElementHessians(_edgeEdgeCollisions.size());
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];
    VECTOR3 e0n = (vs[1] - vs[0]).normalized();
    VECTOR3 e1n = (vs[3] - vs[2]).normalized();
    if (abs(e0n.dot(e1n)) > 1 - 1e-2) {//near-parallel guard
      perElementHessians[i] = MATRIX12::Zero();
      continue;  
    }

    bool isNeg = eeEnergy.getLen(eeEnergy.makeState(vs)) < 0;
    vector<int> extraInfo(1);
    extraInfo[0] = 1;
    if (isNeg != _edgeEdgeIntersections[i]){
      extraInfo[0] = -1;
    }

    //Hessian calculation
    const MATRIX12 H = -_edgeEdgeCollisionAreas[i] * eeEnergy.clampedHessian(vs, extraInfo);
    perElementHessians[i] = H;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const MATRIX12& H = perElementHessians[i];
    const VECTOR2I& edge0 = _edges[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edges[_edgeEdgeCollisions[i].second];

    vector<int> vertexIndices(4);
    vertexIndices[0] = edge0[0];
    vertexIndices[1] = edge0[1];
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];

    for (int y = 0; y < 4; y++)
    {
      int yVertex = vertexIndices[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = vertexIndices[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }
 
  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the collision force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TRIANGLE_MESH::computeVertexFaceCollisionClampedHessian(const ENERGY_12D& vfEnergy) const
{
  TIMER functionTimer(__FUNCTION__);

  vector<MATRIX12> perElementHessians(_vertexFaceCollisionTets.size());
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    vector<VECTOR3> vs(4);
    for (unsigned int j = 0; j < 4; j++)
      vs[j] = _vertices[_vertexFaceCollisionTets[i][j]];
    const MATRIX12 H = -_vertexFaceCollisionAreas[i] * vfEnergy.clampedHessian(vs);
    perElementHessians[i] = H;
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _vertexFaceCollisionTets.size(); i++)
  {
    const VECTOR4I& tet = _vertexFaceCollisionTets[i];
    const MATRIX12& H = perElementHessians[i];
    for (int y = 0; y < 4; y++)
    {
      int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(3 * xVertex + a, 3 * yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// compute whether one vertex is inside the vertex one right of another
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeSurfaceVertexOneRings()
{
  _insideSurfaceVertexOneRing.clear();
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    // gonna be lazy for a moment here
    const VECTOR2I edge = _edges[x];
    _insideSurfaceVertexOneRing[pair<int, int>(edge[0], edge[1])] = true;
    _insideSurfaceVertexOneRing[pair<int, int>(edge[1], edge[0])] = true;
  }
}

///////////////////////////////////////////////////////////////////////
// build a consistent tet/flap ordering from two surface triangles
///////////////////////////////////////////////////////////////////////
VECTOR4I TRIANGLE_MESH::buildFlap(const int triangle0, const int triangle1) const
{
  // they're legal indices, right?
  assert(triangle0 >= 0);
  assert(triangle1 >= 0);
  assert(triangle0 < (int)_triangles.size());
  assert(triangle1 < (int)_triangles.size());

  // they are neighbors, right?
  assert(areTriangleNeighbors(triangle0, triangle1));

  const VECTOR3I f0 = _triangles[triangle0];
  const VECTOR3I f1 = _triangles[triangle1];

  int firstMatch = -1;
  int secondMatch = -1;
  int unmatched0 = -1;
  int unmatched1 = -1;

  // find the tet indices for the first face
  for (int x = 0; x < 3; x++)
  {
    // let's search for this index
    int i0 = f0[x];

    bool matchFound = false;
    for (int y = 0; y < 3; y++)
    {
      // see if it matches
      if (i0 == f1[y])
      {
        // which match is it?
        if (firstMatch == -1)
          firstMatch = i0;
        else
          secondMatch = i0;
        
        matchFound = true;
      }
    }

    if (!matchFound)
      unmatched0 = i0;
  }

  // find the unmatched vertex from the second face
  for (int x = 0; x < 3; x++)
  {
    // let's search for this index
    int i1 = f1[x];
    if (i1 != firstMatch && i1 != secondMatch)
      unmatched1 = i1;
  }

  // we did find one, right?
  assert(unmatched1 != -1);

  // build a tet/flap
  ///////////////////////////
  //                       //                                                                         
  //           1           //
  //                       //
  //           o           //
  //          /|\          //
  //         / | \         //
  //        /  |  \        //
  //    0  o   |   o  3    //
  //        \  |  /        //
  //         \ | /         //
  //          \|/          //
  //           o           //
  //                       //
  //           2           //
  //                       //
  ///////////////////////////
  VECTOR4I tet;
  tet[0] = unmatched0;
  tet[1] = secondMatch;
  tet[2] = firstMatch;
  tet[3] = unmatched1;

  return tet;
}

///////////////////////////////////////////////////////////////////////
// compute flaps and their areas
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeFlaps()
{
  _flaps.clear();
  for (unsigned int x = 0; x < _edgeTriangleNeighbors.size(); x++)
  {
    const VECTOR2I& edgeTriangles = _edgeTriangleNeighbors[x];
    if (edgeTriangles[0] == -1 || edgeTriangles[1] == -1) continue;

    const VECTOR4I flap = buildFlap(edgeTriangles[0], 
                                    edgeTriangles[1]);
    _flaps.push_back(flap);

    const REAL area = _restTriangleAreas[edgeTriangles[0]] + 
                      _restTriangleAreas[edgeTriangles[1]];
    _restFlapAreas.push_back(area);
  }
}

///////////////////////////////////////////////////////////////////////
// get mass-weighted global translation
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::getTranslation() const
{
  VECTOR3 vertexSum;
  vertexSum.setZero();

  REAL areaSum = 0;
  assert(_vertices.size() == _restOneRingAreas.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    areaSum += _restOneRingAreas[x];
    vertexSum += _vertices[x] * _restOneRingAreas[x];
  }

  return vertexSum * (1.0 / areaSum);
}

///////////////////////////////////////////////////////////////////////
// get volume-weighted global translation, for the rest state
///////////////////////////////////////////////////////////////////////
VECTOR3 TRIANGLE_MESH::getRestTranslation() const
{
  VECTOR3 vertexSum;
  vertexSum.setZero();

  REAL areaSum = 0;
  assert(_vertices.size() == _restOneRingAreas.size());
  for (unsigned int x = 0; x < _restVertices.size(); x++)
  {
    areaSum += _restOneRingAreas[x];
    vertexSum += _restVertices[x] * _restOneRingAreas[x];
  }

  return vertexSum * (1.0 / areaSum);
}

///////////////////////////////////////////////////////////////////////
// are these two surface triangles neighbors?
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::areTriangleNeighbors(const int id0, const int id1) const
{
  assert(_triangleNeighbors.size() > 0);
  assert(id0 < (int)_triangleNeighbors.size());
  assert(id1 < (int)_triangleNeighbors.size());

  const VECTOR3I neighbors0 = _triangleNeighbors[id0];

  for (int x = 0; x < 3; x++)
    if (neighbors0[x] == id1)
      return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// compute which vertices are attached to inverted tets
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::computeInvertedVertices()
{
  TIMER functionTimer(__FUNCTION__);
  // first set them all to false
  _invertedVertices.resize(_vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _invertedVertices[x] = false;

  for (unsigned int x = 0; x < _triangles.size(); x++)
  {
    // if the tet is not inverted, move on
    if (invariant3(_Fs[x]) > 0.0)
      continue;

    // if tet is inverted, tags all its vertices
    for (int y = 0; y < 3; y++)
      _invertedVertices[_triangles[x][y]] = true;
  }

  int totalInverted = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    if (_invertedVertices[x])
      totalInverted++;

  //cout << " Total inverted vertices: " << totalInverted << endl;
}

///////////////////////////////////////////////////////////////////////
// set collision pairs (for replays)
///////////////////////////////////////////////////////////////////////
void TRIANGLE_MESH::setCollisionPairs(const vector<pair<int,int> >& vertexFace, 
                                      const vector<pair<int,int> >& edgeEdge)
{
  _vertexFaceCollisions = vertexFace;
  _edgeEdgeCollisions = edgeEdge;
}

///////////////////////////////////////////////////////////////////////
// see if a current surface triangle has been crushed to degeneracy
///////////////////////////////////////////////////////////////////////
bool TRIANGLE_MESH::surfaceTriangleIsDegenerate(const int surfaceTriangleID)
{
  assert(surfaceTriangleID >= 0);
  assert(surfaceTriangleID < (int)_triangles.size());

  // get the rest area
  vector<VECTOR3> vertices(3);
  vertices[0] = _restVertices[_triangles[surfaceTriangleID][0]];
  vertices[1] = _restVertices[_triangles[surfaceTriangleID][1]];
  vertices[2] = _restVertices[_triangles[surfaceTriangleID][2]];
  const REAL restArea = triangleArea(vertices);

  // get the deformed area
  vertices[0] = _vertices[_triangles[surfaceTriangleID][0]];
  vertices[1] = _vertices[_triangles[surfaceTriangleID][1]];
  vertices[2] = _vertices[_triangles[surfaceTriangleID][2]];
  const REAL deformedArea = triangleArea(vertices);

  const REAL relativeArea = deformedArea / restArea;

  const REAL degeneracyEps = 1e-4;
  if (relativeArea < degeneracyEps) return true;

  return false;
}

///////////////////////////////////////////////////////////////////////
// get Procrustes-style global rotation using Eqn. 7 and surrounding
// text from Muller et al's "Meshless Deformations Based on Shape 
// Matching" from SIGGRAPH 2005
///////////////////////////////////////////////////////////////////////
MATRIX3 TRIANGLE_MESH::getRotation() const
{
  // trying to follow Muller's notation here
  const VECTOR3 x_cm0 = getRestTranslation();
  const VECTOR3 x_cm  = getTranslation();

  // left matrix in Eqn. 7 of Muller's paper
  MATRIX3 Apq;
  Apq.setZero();
  for (unsigned int x = 0; x < _restVertices.size(); x++)
  {
    const VECTOR3 p = _vertices[x] - x_cm;
    const VECTOR3 q = _restVertices[x] - x_cm0;
    Apq += _restOneRingAreas[x] * (p * q.transpose());
  }

  // get the rotation
  MATRIX3 R,S;
  polarDecomposition(Apq, R, S);
  return R;
}

}
