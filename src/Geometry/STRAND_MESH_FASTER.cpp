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
#include "STRAND_MESH_FASTER.h"
#include "util/TIMER.h"
#include "LINE_INTERSECT.h"

namespace HOBAK {

using namespace std;

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
STRAND_MESH_FASTER::STRAND_MESH_FASTER(const vector<VECTOR3>& restVertices,
                                       const vector<vector<int> >& strandIndices,
                                       const REAL& E,        // Young's modulus
                                       const REAL& nu,       // Poissons' ratio
                                       const REAL& density,
                                       const REAL& radiusA,
                                       const REAL& radiusB) :
  STRAND_MESH(restVertices, strandIndices, E, nu, density, radiusA, radiusB)
{
  _collisionTree = new AABB_TREE(_vertices, &_edgeIndices);

  _perBendHessians.resize(_totalBends);
  _perEdgeHessians.resize(_totalEdges);

  computeMatrixSparsity();
  computeCompressedIndices();

  //_collisionEps = 0;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
STRAND_MESH_FASTER::~STRAND_MESH_FASTER()
{
  if (_collisionTree)
    delete _collisionTree;
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::computeStretchingHessian(const STRAND::STRETCHING& stretching) const
{
  vector<MATRIX6> perEdgeHessians(_totalEdges);
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    //const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialClampedHessian(positions, 1.0 / _restEdgeLengths[i]);
    const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialHessian(positions, 1.0 / _restEdgeLengths[i]);
    perEdgeHessians[i] = -1.0 * H;
  }
  
  return buildPerEdgeMatrix(perEdgeHessians);
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::computeStretchingClampedHessian(const STRAND::STRETCHING& stretching) const
{
  TIMER functionTimer(string("STRAND_MESH_FASTER::") + string(__FUNCTION__));
  vector<MATRIX6> perEdgeHessians(_totalEdges);

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialClampedHessian(positions, 1.0 / _restEdgeLengths[i]);
    perEdgeHessians[i] = -1.0 * clampEigenvalues(H);
  }

  return buildPerEdgeMatrix(perEdgeHessians);
}

///////////////////////////////////////////////////////////////////////
// bending force gradients
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::computeTanBendingHessian()
{
  TIMER functionTimer(string("STRAND_MESH_FASTER::") + string(__FUNCTION__));
  vector<MATRIX11> perBendHessians(_totalBends);
  if (_bendingForceFilterEnabled)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " BENDING FORCES FILTERED " << endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  }

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const MATRIX2& B = _Bs[x];
    const REAL len = _refVertexLengths[x];

    const VECTOR2& kappa = _kappas[x];
    const VECTOR2& kappaBar = _kappaBars[x];
    const MATRIX11x2& gradKappa = computeGradKappa(x);

    MATRIX11 localJ = 1.0 / len * gradKappa * B * gradKappa.transpose();

    const pair<MATRIX11, MATRIX11>& hessKappa = computeHessianKappa(x);
    const VECTOR2 temp = 1.0 / len * (kappa - kappaBar).transpose() * B;
    localJ += temp(0) * hessKappa.first + temp(1) * hessKappa.second;
    perBendHessians[x] = -1.0 * localJ;

    if (_bendingForceFilterEnabled)
    {
      if (kappa.norm() > _bendingFilterThreshold)
        perBendHessians[x] *= 0.0;
    }
  }

  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision Hessians using x-based formulation
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::computeEdgeEdgeCollisionClampedHessian(const ENERGY_12D& edgeEdgeEnergy) const
{
  TIMER functionTimer(string("STRAND_MESH_FASTER::") + string(__FUNCTION__));

  vector<MATRIX12> perEdgeHessians(_edgeEdgeCollisions.size());
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edgeIndices[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edgeIndices[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];

#if 1
    vector<int> extraInfo;
    extraInfo.push_back(1);
    REAL len = edgeEdgeEnergy.getLen(edgeEdgeEnergy.makeState(vs));
    if(len < 0){
      extraInfo[0] = -1;
    }
    const MATRIX12 H = -_edgeEdgeCollisionAreas[i] * edgeEdgeEnergy.clampedHessian(vs,extraInfo);
#else
    const VECTOR2& a = _edgeEdgeCoordinates[i].first;
    const VECTOR2& b = _edgeEdgeCoordinates[i].second;
    const MATRIX12 H = -_edgeEdgeCollisionAreas[i] * _edgeEdgeEnergy->clampedHessian(vs,a,b);
#endif
    perEdgeHessians[i] = H;
  }

  return buildEdgeEdgeMatrix(perEdgeHessians);
}

#if 0
///////////////////////////////////////////////////////////////////////
// shared matrix construction between clamped and unclamped
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::buildPerBendMatrix(const vector<MATRIX11>& perBendHessians) const
{
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  
  // do the vertex DOFs first
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];
    const MATRIX11& H = perBendHessians[x];

    for (unsigned int i = 0; i < 3; i++)
    {
      const unsigned int vi = vertexIndices[i];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
           {
              const REAL entry = H(4 * i + a, 4 * j + b);
              TRIPLET triplet(3 * vi + a, 3 * vj + b,entry);
              triplets.push_back(triplet);
           }
      }
    }
  }

  // do the edge-edge (twist) DOFs second
  const unsigned int vertexEnd = 3 * _totalVertices;
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];
    const MATRIX11& H = perBendHessians[x];

    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      for (unsigned int j = 0; j < 2; j++)
      {
        const unsigned int ej = edgeIndices[j];
        TRIPLET triplet(vertexEnd + ei, vertexEnd + ej, H(4 * i + 3, 4 * j + 3));
        triplets.push_back(triplet);
      }
    }
  }

  // edge-vertex is last
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];
    const VECTOR2I& edgeIndices = _bendEdges[x];
    const MATRIX11& H = perBendHessians[x];

    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        for (unsigned int a = 0; a < 3; a++)
        {
          const unsigned int row = vertexEnd + ei;
          const unsigned int col = 3 * vj + a;
          const REAL& entry = H(4 * i + 3, 4 * j + a);
          TRIPLET triplet0(row,col,entry);
          TRIPLET triplet1(col,row,entry);
          triplets.push_back(triplet0);
          triplets.push_back(triplet1);
        }
      }
    }
  }

  const int DOFs = 3 * _totalVertices + _totalEdges;
  SPARSE_MATRIX result(DOFs, DOFs);
  result.setFromTriplets(triplets.begin(), triplets.end());

  return result;
}
#endif

///////////////////////////////////////////////////////////////////////
// bending force gradients
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::computeTanBendingClampedHessian()
{
  TIMER functionTimer(string("STRAND_MESH_FASTER::") + string(__FUNCTION__));
  vector<MATRIX11> perBendHessians(_totalBends);
  if (_bendingForceFilterEnabled)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " BENDING FORCES FILTERED " << endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  }

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const MATRIX2& B = _Bs[x];
    const REAL len = _refVertexLengths[x];

    const VECTOR2& kappa = _kappas[x];
    const VECTOR2& kappaBar = _kappaBars[x];
    const MATRIX11x2& gradKappa = computeGradKappa(x);

    MATRIX11 localJ = 1.0 / len * gradKappa * B * gradKappa.transpose();

    const pair<MATRIX11, MATRIX11>& hessKappa = computeHessianKappa(x);
    const VECTOR2 temp = 1.0 / len * (kappa - kappaBar).transpose() * B;
    localJ += temp(0) * hessKappa.first + temp(1) * hessKappa.second;

    perBendHessians[x] = -1.0 * clampEigenvalues(localJ);
    if (_bendingForceFilterEnabled)
    {
      if (kappa.norm() > _bendingFilterThreshold)
        perBendHessians[x] *= 0.0;
    }
  }

  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::computeTanTwistingHessian() const
{
  TIMER functionTimer(string("STRAND_MESH_FASTER::") + string(__FUNCTION__));
  vector<MATRIX11> perBendHessians(_totalBends);

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    const VECTOR11& gradTwist = computeTanGradTwist(x);
    const MATRIX11& hessTwist = computeTanHessianTwist(x);

    MATRIX11 localJ = kt / len * ((twist - undeformedTwist) * hessTwist
                       + gradTwist * gradTwist.transpose());
    perBendHessians[x] = -1.0 * localJ;
  }

  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::computeTanTwistingClampedHessian() const
{
  TIMER functionTimer(string("STRAND_MESH_FASTER::") + string(__FUNCTION__));
  vector<MATRIX11> perBendHessians(_totalBends);

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    const VECTOR11& gradTwist = computeTanGradTwist(x);
    const MATRIX11& hessTwist = computeTanHessianTwist(x);

    MATRIX11 localJ = kt / len * ((twist - undeformedTwist) * hessTwist
                       + gradTwist * gradTwist.transpose());
    perBendHessians[x] = -1.0 * clampEigenvalues(localJ);
  }

  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
// edge-edge collision detection
//
// TODO: each edge only owns one vertex; the downstream one
///////////////////////////////////////////////////////////////////////
void STRAND_MESH_FASTER::computeEdgeEdgeCollisions(const REAL& collisionEps, const bool verbose)
{
#if 0
  STRAND_MESH::computeEdgeEdgeCollisions();
  return;
#endif

  TIMER functionTimer(string("STRAND_MESH_FASTER::") + string(__FUNCTION__));
  //const REAL collisionEps = edgeEdgeEnergy.eps();

  // for visualization purposes
  _edgeEdgeCollisionsOld = _edgeEdgeCollisions;
  _edgeEdgeCoordinatesOld = _edgeEdgeCoordinates;

  _edgeEdgeCollisions.clear();
  _edgeEdgeIntersections.clear();
  _edgeEdgeCoordinates.clear();
  _edgeEdgeCollisionAreas.clear();

  _collisionTree->refit();
  // get the nearest edge to each edge, not including itself
  // and ones where it shares a vertex
  for (unsigned int x = 0; x < _edgeIndices.size(); x++)
  {
    //int closestEdge = -1;
    //REAL closestDistance = FLT_MAX;
    VECTOR2 aClosest(-1,-1);
    VECTOR2 bClosest(-1,-1);
    const VECTOR2I& outerEdge = _edgeIndices[x];
    const VECTOR3& v0 = _vertices[outerEdge[0]];
    const VECTOR3& v1 = _vertices[outerEdge[1]];
    //const unsigned int outerFlat = outerEdge[0] + outerEdge[1] * _edgeIndices.size();

    vector<int> nearbyEdges;
    //_collisionTree->nearbyEdges(_edgeIndices[x], _collisionEps, nearbyEdges);
    _collisionTree->nearbyEdges(_edgeIndices[x], collisionEps, nearbyEdges);

    // find the closest other edge
    for (unsigned int y = 0; y < nearbyEdges.size(); y++)
    {
      if ((int)x == nearbyEdges[y]) continue;

      // skip if index is smaller -- don't want to double count nearby edges
      // (a,b) and (b,a)
      if ((unsigned int)nearbyEdges[y] < x) continue;

      const VECTOR2I innerEdge = _edgeIndices[nearbyEdges[y]];

      // if they share a vertex, skip it
      if ((outerEdge[0] == innerEdge[0]) || (outerEdge[0] == innerEdge[1]) ||
          (outerEdge[1] == innerEdge[0]) || (outerEdge[1] == innerEdge[1]))
        continue;

      const VECTOR3& v2 = _vertices[innerEdge[0]];
      const VECTOR3& v3 = _vertices[innerEdge[1]];

      VECTOR3 innerPoint, outerPoint;
      IntersectLineSegments(v0, v1, v2, v3,
                            outerPoint, innerPoint);

      const REAL distance = (innerPoint - outerPoint).norm();

      // if it's not close enough, skip it, but if it is close enough,
      // it's fine to add multiple contacts
      //if (distance > _collisionEps) continue;
      if (distance > collisionEps) continue;

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
      //const REAL skipEps = 1e-4;
      const REAL skipEps = 0;
      if ((a[0] < skipEps) || (a[0] > 1.0 - skipEps)) continue;
      if ((a[1] < skipEps) || (a[1] > 1.0 - skipEps)) continue;
      if ((b[0] < skipEps) || (b[0] > 1.0 - skipEps)) continue;
      if ((b[1] < skipEps) || (b[1] > 1.0 - skipEps)) continue;

      // sanity check, what's the difference found here?
      const VECTOR3 middle0 = a[0] * _vertices[outerEdge[0]] + a[1] * _vertices[outerEdge[1]];
      const VECTOR3 middle1 = b[0] * _vertices[innerEdge[0]] + b[1] * _vertices[innerEdge[1]];
      const VECTOR3 diff = middle0 - middle1;

      //if (diff.norm() > _collisionEps) continue;
      if (diff.norm() > collisionEps) continue;
      //cout << " collision diff: " << diff.norm() << endl;

      pair<int,int> collision(x, nearbyEdges[y]);
      _edgeEdgeCollisions.push_back(collision);

      pair<VECTOR2,VECTOR2> coordinate(a, b);
      _edgeEdgeCoordinates.push_back(coordinate);

      // get the areas too
      //const VECTOR2I innerEdge = _edgeIndices[nearbyEdges[y]];
      const pair<int,int> outerPair(outerEdge[0], outerEdge[1]);
      const pair<int,int> innerPair(innerEdge[0], innerEdge[1]);
      const REAL xArea = _restEdgeLengths[_edgeHash[outerPair]];
      const REAL closestArea = _restEdgeLengths[_edgeHash[innerPair]];
      _edgeEdgeCollisionAreas.push_back(xArea + closestArea);
    }
  }
  assert(_edgeEdgeCollisions.size() == _edgeEdgeCoordinates.size());
  //if (_edgeEdgeCollisions.size() > 0)
  if (verbose)
    cout << " Found " << _edgeEdgeCollisions.size() << " edge-edge collisions " << endl;

//#define VERY_VERBOSE 1
#if 0
  for (unsigned int x = 0; x < _edgeEdgeCollisions.size(); x++)
  {
    const pair<int,int> collision = _edgeEdgeCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")";

    const VECTOR2& coordinate0 = _edgeEdgeCoordinates[x].first;
    const VECTOR2& coordinate1 = _edgeEdgeCoordinates[x].second;
    cout << "[ " << coordinate0.transpose() << "] [" << coordinate1.transpose() << "]" << endl;
  }
#endif
}

///////////////////////////////////////////////////////////////////////
// fast computation of all the clamped Hessians for the elasticity 
// energies
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH_FASTER::computeClampedElasticityHessian(const STRAND::STRETCHING& stretching)
{
  TIMER computeTimer("computeClampedElasticityHessian, compute");

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    // wipe the previous results
    _perBendHessians[x].setZero();

    // bending energy first
    const MATRIX2& B = _Bs[x];
    const REAL len = _refVertexLengths[x];

    const VECTOR2& kappa = _kappas[x];
    const VECTOR2& kappaBar = _kappaBars[x];
    const MATRIX11x2& gradKappa = computeGradKappa(x);

    MATRIX11 bendJ = 1.0 / len * gradKappa * B * gradKappa.transpose();

    const pair<MATRIX11, MATRIX11>& hessKappa = computeHessianKappa(x);
    const VECTOR2 temp = 1.0 / len * (kappa - kappaBar).transpose() * B;
    bendJ += temp(0) * hessKappa.first + temp(1) * hessKappa.second;

    _perBendHessians[x] = -1.0 * clampEigenvalues(bendJ);
    if (_bendingForceFilterEnabled)
    {
      if (kappa.norm() > _bendingFilterThreshold)
        _perBendHessians[x] *= 0.0;
    }
    
    // then the twisting energy
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    const VECTOR11& gradTwist = computeTanGradTwist(x);
    const MATRIX11& hessTwist = computeTanHessianTwist(x);

    MATRIX11 twistJ = kt / len * ((twist - undeformedTwist) * hessTwist
                       + gradTwist * gradTwist.transpose());
    _perBendHessians[x] += -1.0 * clampEigenvalues(twistJ);
  }

  // compute stretching Hessian
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialClampedHessian(positions, 1.0 / _restEdgeLengths[i]);
    _perEdgeHessians[i] = -1.0 * clampEigenvalues(H);
  }
  computeTimer.stop();

  // push everything into _sparseK
  TIMER assemblyTimer("computeClampedElasticityHessian, assembly");
  const unsigned int nonZeros = _sparseK.nonZeros();
  REAL* base = _sparseK.valuePtr();
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < nonZeros; x++)
  {
    // zeroing out here instead of using a global setZero, because that will
    // erase the sparsity pattern
    base[x] = 0;

    // do the per-bend, bending and twisting first
    const vector<VECTOR3I>& bendGather = _hessianBendGathers[x];
    for (unsigned int y = 0; y < bendGather.size(); y++)
    {
      const VECTOR3I& lookup = bendGather[y];
      const int& bendIndex = lookup[0];
      const int& row = lookup[1];
      const int& col = lookup[2];

      assert(bendIndex < (int)_perBendHessians.size());

      assert(row < 11);
      assert(col < 11);

      base[x] += _perBendHessians[bendIndex](row, col);
    }
    
    // then the per-edge stretching
    const vector<VECTOR3I>& edgeGather = _hessianEdgeGathers[x];
    for (unsigned int y = 0; y < edgeGather.size(); y++)
    {
      const VECTOR3I& lookup = edgeGather[y];
      const int& edgeIndex = lookup[0];
      const int& row = lookup[1];
      const int& col = lookup[2];

      assert(edgeIndex < (int)_perEdgeHessians.size());

      assert(row < 11);
      assert(col < 11);

      base[x] += _perEdgeHessians[edgeIndex](row, col);
    }
  }
  assemblyTimer.stop();

  //return _sparseK;
  return _sparseK.selfadjointView<Eigen::Lower>();
}

///////////////////////////////////////////////////////////////////////
// bake out the sparsity pattern
///////////////////////////////////////////////////////////////////////
void STRAND_MESH_FASTER::computeMatrixSparsity()
{
  if (_edgeEnd)
  {
    computeMatrixSparsityEdgeEnd();
    return;
  }
  computeMatrixSparsityInterleaved();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH_FASTER::computeMatrixSparsityEdgeEnd()
{
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  
  // do the vertex DOFs first
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];
    for (unsigned int i = 0; i < 3; i++)
    {
      const unsigned int vi = vertexIndices[i];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
           {
              // just set it to one
              TRIPLET triplet(3 * vi + a, 3 * vj + b, 1);
              triplets.push_back(triplet);
           }
      }
    }
  }

  // do the edge-edge (twist) DOFs second
  const unsigned int vertexEnd = 3 * _totalVertices;
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];
    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      for (unsigned int j = 0; j < 2; j++)
      {
        // just set it to one
        const unsigned int ej = edgeIndices[j];
        TRIPLET triplet(vertexEnd + ei, vertexEnd + ej, 1.0);
        triplets.push_back(triplet);
      }
    }
  }

  // edge-vertex is last
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];
    const VECTOR2I& edgeIndices = _bendEdges[x];

    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        for (unsigned int a = 0; a < 3; a++)
        {
          const unsigned int row = vertexEnd + ei;
          const unsigned int col = 3 * vj + a;
          TRIPLET triplet0(row,col,1);
          TRIPLET triplet1(col,row,1);
          triplets.push_back(triplet0);
          triplets.push_back(triplet1);
        }
      }
    }
  }

  const int DOFs = 3 * _totalVertices + _totalEdges;
  _sparseK = SPARSE_MATRIX(DOFs, DOFs);
  _sparseK.setFromTriplets(triplets.begin(), triplets.end());
  _sparseK.makeCompressed();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH_FASTER::computeMatrixSparsityInterleaved()
{
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  
  // do the vertex DOFs first
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];
    for (unsigned int i = 0; i < 3; i++)
    {
      const unsigned int vi = vertexIndices[i];
      const unsigned int viGlobal = _globalVertexIndices[vi];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        const unsigned int vjGlobal = _globalVertexIndices[vj];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
           {
              // only do lower triangle
              const int row = viGlobal + a;
              const int col = vjGlobal + b;
              if (col > row) continue;
              TRIPLET triplet(row, col, 1);
              triplets.push_back(triplet);

              // just set it to one
              //TRIPLET triplet(viGlobal + a, vjGlobal + b, 1);
              //triplets.push_back(triplet);
           }
      }
    }
  }

  // do the edge-edge (twist) DOFs second
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];
    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      const unsigned int eiGlobal = _globalEdgeIndices[ei];
      for (unsigned int j = 0; j < 2; j++)
      {
        // just set it to one
        const unsigned int ej = edgeIndices[j];
        const unsigned int ejGlobal = _globalEdgeIndices[ej];
        
        // only do lower triangle
        const int row = eiGlobal;
        const int col = ejGlobal;
        if (col > row) continue;
        
        TRIPLET triplet(eiGlobal, ejGlobal, 1.0);
        triplets.push_back(triplet);
      }
    }
  }

  // edge-vertex is last
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];
    const VECTOR2I& edgeIndices = _bendEdges[x];

    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      const unsigned int eiGlobal = _globalEdgeIndices[ei];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        const unsigned int vjGlobal = _globalVertexIndices[vj];
        for (unsigned int a = 0; a < 3; a++)
        {
          const unsigned int row = eiGlobal;
          const unsigned int col = vjGlobal + a;
          TRIPLET triplet0(row,col,1);
          TRIPLET triplet1(col,row,1);
          // only do lower triangle
          if (row > col)
            triplets.push_back(triplet0);
          // only do lower triangle
          if (col > row)
            triplets.push_back(triplet1);
        }
      }
    }
  }

  const int DOFs = 3 * _totalVertices + _totalEdges;
  _sparseK = SPARSE_MATRIX(DOFs, DOFs);
  _sparseK.setFromTriplets(triplets.begin(), triplets.end());
  _sparseK.makeCompressed();
}

///////////////////////////////////////////////////////////////////////
// find the compressed index mapping
///////////////////////////////////////////////////////////////////////
void STRAND_MESH_FASTER::computeCompressedIndices()
{
  if (_edgeEnd)
  {
    computeCompressedIndicesEdgeEnd();
    return;
  }
  computeCompressedIndicesInterleaved();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH_FASTER::computeCompressedIndicesEdgeEnd()
{
  TIMER functionTimer(__FUNCTION__);

  cout << " Hashing indices ... " << flush;

  // cache the beginning of the storage
  REAL* base = _sparseK.valuePtr();

  for (unsigned int x = 0; x < _sparseK.outerSize(); x++)
  {
    for (SPARSE_MATRIX::InnerIterator it(_sparseK, x); it; ++it)
    {
      // make the (row, col) pair
      const pair<int, int> rowCol(it.row(), it.col());

      // get the index
      const int index = (int)(&it.value() - base);

      // find the address and store it in the map
      _compressedIndex[rowCol] = index;
    }
  }
  cout << "done." << endl;

  cout << " Computing compressed bend indices ... " << flush;
  // allocate an array for each non-zero matrix entry
  _hessianBendGathers.resize(_sparseK.nonZeros());

  // do the vertex DOFs first
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];

    for (unsigned int i = 0; i < 3; i++)
    {
      const unsigned int vi = vertexIndices[i];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            // do the lookup, see where this is stored globally
            const pair<int, int> rowCol(3 * vi + a, 3 * vj + b);
            const auto iter = _compressedIndex.find(rowCol);
            const int index = iter->second;

            // store the entry and H this corresponds to
            VECTOR3I vertexMapping;
            vertexMapping[0] = x;
            vertexMapping[1] = 4 * i + a;
            vertexMapping[2] = 4 * j + b;
            _hessianBendGathers[index].push_back(vertexMapping);
          }
      }
    }
  }

  // do the edge-edge (twist) DOFs second
  const unsigned int vertexEnd = 3 * _totalVertices;
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];
    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      for (unsigned int j = 0; j < 2; j++)
      {
        const unsigned int ej = edgeIndices[j];
        // do the lookup, see where this is stored globally
        const pair<int, int> rowCol(vertexEnd + ei, vertexEnd + ej);
        const auto iter = _compressedIndex.find(rowCol);
        const int index = iter->second;

        // store the entry and H this corresponds to
        VECTOR3I edgeMapping;
        edgeMapping[0] = x;
        edgeMapping[1] = 4 * i + 3;
        edgeMapping[2] = 4 * j + 3;
        _hessianBendGathers[index].push_back(edgeMapping);
      }
    }
  }

  // edge-vertex is last
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];
    const VECTOR2I& edgeIndices = _bendEdges[x];
    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        for (unsigned int a = 0; a < 3; a++)
        {
          const unsigned int row = vertexEnd + ei;
          const unsigned int col = 3 * vj + a;

          // do the lookup, see where this is stored globally
          const pair<int, int> rowCol(row, col);
          const auto iterRowCol = _compressedIndex.find(rowCol);
          const int indexRowCol = iterRowCol->second;
          
          const pair<int, int> colRow(col, row);
          const auto iterColRow= _compressedIndex.find(colRow);
          const int indexColRow = iterColRow->second;

          // store the entry and H this corresponds to
          VECTOR3I rowColMapping;
          rowColMapping[0] = x;
          rowColMapping[1] = 4 * i + 3;
          rowColMapping[2] = 4 * j + a;
          _hessianBendGathers[indexRowCol].push_back(rowColMapping);

          VECTOR3I colRowMapping;
          colRowMapping[0] = x;
          colRowMapping[1] = 4 * i + 3;
          colRowMapping[2] = 4 * j + a;
          _hessianBendGathers[indexColRow].push_back(colRowMapping);
        }
      }
    }
  }
  cout << "done." << endl;

  cout << " Computing compressed edge indices ... " << flush;

  // allocate an array for each non-zero matrix entry
  _hessianEdgeGathers.resize(_sparseK.nonZeros());
  
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];
    for (int y = 0; y < 2; y++)
    {
      const int yVertex = edge[y];
      for (int x = 0; x < 2; x++)
      {
        const int xVertex = edge[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const unsigned int row = 3 * xVertex + a;
            const unsigned int col = 3 * yVertex + b;
          
            // do the lookup, see where this is stored globally
            const pair<int, int> rowCol(row, col);
            const auto iter = _compressedIndex.find(rowCol);
            const int index = iter->second;

            // store the entry and H this corresponds to
            VECTOR3I rowColMapping;
            rowColMapping[0] = i;
            rowColMapping[1] = 3 * x + a;
            rowColMapping[2] = 3 * y + b;
            _hessianEdgeGathers[index].push_back(rowColMapping);
          }
      }
    }
  }
  cout << "done." << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH_FASTER::computeCompressedIndicesInterleaved()
{
  TIMER functionTimer(__FUNCTION__);

  cout << " Hashing indices ... " << flush;

  // cache the beginning of the storage
  REAL* base = _sparseK.valuePtr();

  for (unsigned int x = 0; x < _sparseK.outerSize(); x++)
  {
    for (SPARSE_MATRIX::InnerIterator it(_sparseK, x); it; ++it)
    {
      // only do lower triangle
      if (it.col() > it.row()) continue;

      // make the (row, col) pair
      const pair<int, int> rowCol(it.row(), it.col());

      // get the index
      const int index = (int)(&it.value() - base);

      // find the address and store it in the map
      _compressedIndex[rowCol] = index;
    }
  }
  cout << "done." << endl;

  cout << " Computing compressed bend indices ... " << flush;
  // allocate an array for each non-zero matrix entry
  _hessianBendGathers.resize(_sparseK.nonZeros());

  // do the vertex DOFs first
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];

    for (unsigned int i = 0; i < 3; i++)
    {
      const unsigned int vi = vertexIndices[i];
      const unsigned int viGlobal = _globalVertexIndices[vi];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        const unsigned int vjGlobal = _globalVertexIndices[vj];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            // only do lower triangle
            const int row = viGlobal + a;
            const int col = vjGlobal + b;
            if (col > row) continue;
            const pair<int, int> rowCol(row, col);

            // do the lookup, see where this is stored globally
            //const pair<int, int> rowCol(viGlobal + a, vjGlobal + b);
            const auto iter = _compressedIndex.find(rowCol);
            const int index = iter->second;

            // store the entry and H this corresponds to
            VECTOR3I vertexMapping;
            vertexMapping[0] = x;
            vertexMapping[1] = 4 * i + a;
            vertexMapping[2] = 4 * j + b;
            _hessianBendGathers[index].push_back(vertexMapping);
          }
      }
    }
  }

  // do the edge-edge (twist) DOFs second
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];
    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      const unsigned int eiGlobal = _globalEdgeIndices[ei];
      for (unsigned int j = 0; j < 2; j++)
      {
        const unsigned int ej = edgeIndices[j];
        const unsigned int ejGlobal = _globalEdgeIndices[ej];
        
        // only do lower triangle
        const int row = eiGlobal;
        const int col = ejGlobal;
        if (col > row) continue;

        // do the lookup, see where this is stored globally
        const pair<int, int> rowCol(eiGlobal, ejGlobal);
        const auto iter = _compressedIndex.find(rowCol);
        const int index = iter->second;

        // store the entry and H this corresponds to
        VECTOR3I edgeMapping;
        edgeMapping[0] = x;
        edgeMapping[1] = 4 * i + 3;
        edgeMapping[2] = 4 * j + 3;
        _hessianBendGathers[index].push_back(edgeMapping);
      }
    }
  }

  // edge-vertex is last
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& vertexIndices = _bendVertices[x];
    const VECTOR2I& edgeIndices = _bendEdges[x];
    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      const unsigned int eiGlobal = _globalEdgeIndices[ei];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        const unsigned int vjGlobal = _globalVertexIndices[vj];
        for (unsigned int a = 0; a < 3; a++)
        {
          const unsigned int row = eiGlobal;
          const unsigned int col = vjGlobal + a;

          // do the lookup, see where this is stored globally
          const pair<int, int> rowCol(row, col);
          const auto iterRowCol = _compressedIndex.find(rowCol);
          const int indexRowCol = iterRowCol->second;
          
          const pair<int, int> colRow(col, row);
          const auto iterColRow= _compressedIndex.find(colRow);
          const int indexColRow = iterColRow->second;

          // store the entry and H this corresponds to
          VECTOR3I rowColMapping;
          rowColMapping[0] = x;
          rowColMapping[1] = 4 * i + 3;
          rowColMapping[2] = 4 * j + a;

          // only do lower triangle
          if (row > col)
            _hessianBendGathers[indexRowCol].push_back(rowColMapping);

          VECTOR3I colRowMapping;
          colRowMapping[0] = x;
          colRowMapping[1] = 4 * i + 3;
          colRowMapping[2] = 4 * j + a;
          // only do lower triangle
          if (col > row)
            _hessianBendGathers[indexColRow].push_back(colRowMapping);
        }
      }
    }
  }
  cout << "done." << endl;

  cout << " Computing compressed edge indices ... " << flush;

  // allocate an array for each non-zero matrix entry
  _hessianEdgeGathers.resize(_sparseK.nonZeros());
  
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];
    for (int y = 0; y < 2; y++)
    {
      const int yVertex = edge[y];
      const int yIndex = _globalVertexIndices[yVertex];
      for (int x = 0; x < 2; x++)
      {
        const int xVertex = edge[x];
        const int xIndex = _globalVertexIndices[xVertex];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const unsigned int row = xIndex + a;
            const unsigned int col = yIndex + b;
  
            // only do lower triangle
            if (col > row) continue;

            // do the lookup, see where this is stored globally
            const pair<int, int> rowCol(row, col);
            const auto iter = _compressedIndex.find(rowCol);
            const int index = iter->second;

            // store the entry and H this corresponds to
            VECTOR3I rowColMapping;
            rowColMapping[0] = i;
            rowColMapping[1] = 3 * x + a;
            rowColMapping[2] = 3 * y + b;
            _hessianEdgeGathers[index].push_back(rowColMapping);
          }
      }
    }
  }
  cout << "done." << endl;
}

}
