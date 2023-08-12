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
#include "STRAND_NET_MESH.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_E_BENDING.h"
#include "Hyperelastic/Strand/TAN_BENDING.h"
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


namespace HOBAK {
// typedef STRAND::QUADRATIC_E_BENDING BENDING_ENERGY;
typedef STRAND::TAN_BENDING BENDING_ENERGY;

using namespace std;
#define USING_BRUTE_FORCE_CLAMP_STRAND 0
///////////////////////////////////////////////////////////////////////
// accepts a vector of individual strands
///////////////////////////////////////////////////////////////////////
STRAND_NET_MESH::STRAND_NET_MESH(const vector<VECTOR3>& restVertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB) :
    STRAND_MESH()
{
  cout<<"in net initializer.\n";
  _restVertices = restVertices;
  _vertices = restVertices;
  _strandIndices = strandIndices;
  _E = E;
  _G = G;
  _density = density;
  _radiusA = radiusA;
  _radiusB = radiusB;

  // STRAND_MESH
  _totalStrands = _strandIndices.size();

  initialize();
  updateProperties();

  
  cout << " total strands: " << _totalStrands << endl;
  cout << " DOFs:          " << DOFs() << endl;

  // STRAND_MESH_FASTER
  // _collisionTree = new AABB_TREE(_vertices, &_edgeIndices);

  // _perBendHessians.resize(_totalBends);
  // _perEdgeHessians.resize(_totalEdges);

  // computeMatrixSparsity();
  // computeCompressedIndices();
}

STRAND_NET_MESH::STRAND_NET_MESH(const vector<VECTOR3>& restVertices,
              const vector<VECTOR3>& vertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB) :
    STRAND_MESH()
{
  cout<<"in net initializer.\n";
  _vertices = vertices;
  _restVertices = restVertices;
  _strandIndices = strandIndices;
  _E = E;
  _G = G;
  _density = density;
  _radiusA = radiusA;
  _radiusB = radiusB;

  // STRAND_MESH
  _totalStrands = _strandIndices.size();

  initialize();
  updateProperties();
  
  cout << " total strands: " << _totalStrands << endl;
  cout << " DOFs:          " << DOFs() << endl;

  // STRAND_MESH_FASTER
  // _collisionTree = new AABB_TREE(_vertices, &_edgeIndices);

  // _perBendHessians.resize(_totalBends);
  // _perEdgeHessians.resize(_totalEdges);

  // computeMatrixSparsity();
  // computeCompressedIndices();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
STRAND_NET_MESH::~STRAND_NET_MESH()
{
}

REAL STRAND_NET_MESH::clamp(const REAL& input, const REAL& bottom, const REAL& top)
{
  REAL result = (input > bottom) ? input : bottom;
  return result < top ? result : top;
}

VECTOR3 STRAND_NET_MESH::findOrthogonal(const VECTOR3& u)
{
  assert(u.norm() != 0);

  VECTOR3 v;
  v.setZero();
  int max = 0;
  for (int i = 0; i < u.size(); ++i) {
    if (u[i] == 0) {
      v[i] = 1;
      return v;
    }
    if (fabs(u[i]) > fabs(u[max])) max = i;
  }

  int idx = (max + 1) % u.size();
  v[idx] = u[max];
  v[max] = -u[idx];
  v.normalize();

  assert(fabs(u.dot(v)) < 1e-10);

  return v;
}

void STRAND_NET_MESH::initialize()
{
  _vertexToEdge.resize(_restVertices.size());
  for(unsigned int x = 0; x < _restVertices.size(); x++){
    _vertexToEdge.clear();
  }
  for (unsigned int y = 0; y < _totalStrands; y++)
  {
    const vector<int>& strand = _strandIndices[y];

    // need to store this up front, so we can index correctly
    // when making the _bendEdges
    const int edgeStart = _edgeIndices.size();

    // inter-strand bends
    for (unsigned int x = 0; x < strand.size() - 1; x++){
      // first vertex
      for (unsigned int y = 0; y < _vertexToEdge[strand[x]].size(); y++){
        // create the per-bend vertex lists
        const VECTOR2I& refEdge = _vertexToEdge[strand[x]][y];
        const VECTOR3I vertexIndices(refEdge(0), strand[x], strand[x+1]);
        _bendVertices.push_back(vertexIndices);

        // create the per-bend edge lists
        const int edgeOffset = edgeStart + x;
        const VECTOR2I edgeIndices(refEdge(1), edgeOffset);
        _bendEdges.push_back(edgeIndices);
      }

      // second vertex
      for (unsigned int y = 0; y < _vertexToEdge[strand[x + 1]].size(); y++){
        // create the per-bend vertex lists
        const VECTOR2I& refEdge = _vertexToEdge[strand[x + 1]][y];
        const VECTOR3I vertexIndices(refEdge(0), strand[x + 1], strand[x]);
        _bendVertices.push_back(vertexIndices);

        // create the per-bend edge lists
        const int edgeOffset = edgeStart + x;
        const VECTOR2I edgeIndices(refEdge(1), edgeOffset);
        _bendEdges.push_back(edgeIndices);
      }
    }

    // create the edges
    for (unsigned int x = 1; x < strand.size(); x++){
      _edgeIndices.push_back(VECTOR2I(strand[x - 1], strand[x]));
      _vertexToEdge[strand[x - 1]].push_back(VECTOR2I(strand[x], _edgeIndices.size() - 1));
      _vertexToEdge[strand[x]].push_back(VECTOR2I(strand[x - 1], _edgeIndices.size() - 1));
    }

    for (unsigned int x = 0; x < strand.size() - 2; x++)
    {
      // create the per-bend vertex lists
      const VECTOR3I vertexIndices(strand[x], strand[x+1], strand[x+2]);
      _bendVertices.push_back(vertexIndices);

      // create the per-bend edge lists
      const int edgeOffset = edgeStart + x;
      const VECTOR2I edgeIndices(edgeOffset, edgeOffset + 1);
      _bendEdges.push_back(edgeIndices);
    }
  }

  const MATRIX3 I3 = MATRIX3::Identity(), z33 = MATRIX3::Zero();
  const VECTOR3 z3 = VECTOR3::Zero();
  _pFpx << -I3, I3, z33,
           z33, -I3, I3;
  _totalEdges = _edgeIndices.size();
  _totalBends = _bendVertices.size();
  _totalVertices = _restVertices.size();
  // _DOFs = 3 * _totalVertices + _totalEdges;
  _DOFs = 3 * _totalVertices;
  _restThetas.resize(_totalBends);
  // compute rest thetas
  for(unsigned int x = 0; x < _totalBends; x++){
    const VECTOR3I& bendVert = _bendVertices[x];
    const VECTOR3 e0 = _restVertices[bendVert(1)] - _restVertices[bendVert(0)];
    const VECTOR3 e1 = _restVertices[bendVert(2)] - _restVertices[bendVert(1)];
    _restThetas[x] = acos(clamp(e0.dot(e1) / (e0.norm() * e1.norm()), -1.0, 1.0));
  }

  _tangents.resize(_totalEdges);
  _director1s.resize(_totalEdges);
  _director2s.resize(_totalEdges);
  _kbs.resize(_totalBends);
  _edgeLengths.resize(_totalEdges);
  _voronoiLengths.resize(_totalVertices);
  _material1s.resize(_totalEdges);
  _material2s.resize(_totalEdges);

  _thetas.resize(_totalEdges);
  _thetas.setZero();

  _kappas.resize(_totalBends);
  _Bs.resize(_totalBends);
  _refVertexLengths.resize(_totalBends);

  _referenceTwists.resize(_totalBends);
  _twists.resize(_totalBends);
  _vertexMasses.resize(_totalVertices);

  _referenceTwists.setZero();
  _twists.setZero();
  _vertexMasses.setZero();
// cout <<"flag----------------\n";
  computeEdges();
  // cout <<"flag----------------\n";
  computeTangents();
  computeCurvatureBinormals();
  computeEdgeLengths();
  _restEdgeLengths = _edgeLengths;
  computeVoronoiLengths();
// cout <<"flag----------------\n";
  // do the frame walk
  // _director1s[0] = findOrthogonal(_tangents[0]);
  // computeSpaceParallel();

  // computeMaterialDirectors();
  computeVertexMasses();

  // computeKappas();
  // computeTwists();
  computeBs(_E, _radiusA, _radiusB);
  computeRefVertexLengths();

  _kappaBars = _kappas;
  _undeformedTwists = _twists;

  // compute and set the twisting constant
  _kts.resize(_totalBends);

  // kt = (1/8) (E / (1 + nu)) * (a b (a^2 + b^2)
  // or in terms of shear modulus G
  //   kt = (1/4) G * (a b (a^2 + b^2)
  //   G  = (1/2) (E / (1 + nu))
  //const REAL kt = ((1.0 / 8.0) * _E / (1.0 + _nu)) * _radiusA * _radiusB *
  //                (_radiusA * _radiusA + _radiusB * _radiusB);
  //
  // no, this formula for Poisson's ratio is giving weird results. Let's just be more
  // direct and use the shear modulus
  const REAL kt = (M_PI / 4.0) * _G * _radiusA * _radiusB *
                  (_radiusA * _radiusA + _radiusB * _radiusB);
  _kts.setOnes();
  _kts *= kt;
#if 0
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  std::cout << " ZEROED TWISTING" << std::endl;
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  _kts *= 0.0;
#endif

  //_collisionEps = 0.01;
  //_collisionEps = 0.1;
  //_collisionEps = 0;
  //_collisionEps = 1.0;

  // mapping from edge index pairs to _edgeIndices
  for (unsigned int x = 0; x < _edgeIndices.size(); x++)
  {
    pair<int,int> edge(_edgeIndices[x][0], _edgeIndices[x][1]);
    _edgeHash[edge] = x;
  }

  // this gets overwritten by TIMESTEPPER every step, so a dummy is fine
  //const REAL stiffness = 1000.0;
  //_edgeEdgeEnergy = new VOLUME::EDGE_SQRT_COLLISION(stiffness, _collisionEps); // default

  // cache for visualization
  _verticesOld = _vertices;

  // filter bending forces by default, since not doing so causes things to explode
  _bendingForceFilterEnabled = false;
  //_bendingFilterThreshold = 100.0;
  _bendingFilterThreshold = 0.0;
  //_bendingFilterThreshold = 1.0;
  
  // are we putting the edge information at the end?
  _edgeEnd = true;
  // _edgeEnd = false;

  /*
  if (_edgeEnd)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    cout << " Edges are PACKED ON THE END " << endl;
  }
  else
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    cout << " Edges are INTERLEAVED" << endl;
  }
  */

  // build out the indexing into the global vector
  cout<<"in net initialize.\n";
  buildGlobalIndices();
}

void STRAND_NET_MESH::buildGlobalIndices()
{
  _globalVertexIndices.resize(_totalVertices);
  _globalEdgeIndices.resize(0);
  _globalStrandEnds.resize(0);
  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    _globalVertexIndices[x] = 3 * x;
  }
}

void STRAND_NET_MESH::updateProperties()
{
  computeEdges();
  computeTangents();
  computeCurvatureBinormals();
  computeEdgeLengths();
  computeVoronoiLengths();

  // Updating so that density update is taken into account
  computeVertexMasses();

  // computeKappas();
  // computeTwists();
}

const VECTOR STRAND_NET_MESH::getDisplacement() const
{
  VECTOR displacements(_DOFs);

  int index = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++, index++)
      displacements[index] = _vertices[x][y] - _restVertices[x][y];
  }

  return displacements;
}

void STRAND_NET_MESH::setDisplacement(const VECTOR& delta)
{
  assert((unsigned int)delta.size() == _DOFs);

  // back up old positions, for visualization
  _verticesOld = _vertices;

  unsigned int index = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++, index++)
      _vertices[x][y] = _restVertices[x][y] + delta[index];
  }

  updateProperties();
}

VECTOR STRAND_NET_MESH::buildPerBendVector(const vector<VECTOR9>& perBendForces) const
{
  assert(perBendForces.size() == _totalBends);
  VECTOR forces(_DOFs);
  forces.setZero();

  for (unsigned int y = 0; y < _totalBends; y++){
    const VECTOR3I& bendVert = _bendVertices[y];
    const VECTOR9& perBendForce = perBendForces[y];
    for(unsigned int x = 0; x < 3; x++){
      const int vertexID = bendVert[x];
      for(unsigned int i = 0; i < 3; i++){
        forces[vertexID * 3 + i] += perBendForce[x * 3 + i];
      }
    }
  }
  return forces;
}

SPARSE_MATRIX STRAND_NET_MESH::buildPerBendMatrix(const vector<MATRIX9>& perBendHessians) const
{
  assert(perBendHessians.size() == _totalBends);

  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& bendVert = _bendVertices[x];
    const MATRIX9& H = perBendHessians[x];

    for (unsigned int i = 0; i < 3; i++)
    {
      const unsigned int vi = bendVert[i];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = bendVert[j];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
           {
              const REAL entry = H(3 * i + a, 3 * j + b);
              TRIPLET triplet(3 * vi + a, 3 * vj + b,entry);
              triplets.push_back(triplet);
           }
      }
    }
  }

  SPARSE_MATRIX result(_DOFs, _DOFs);
  result.setFromTriplets(triplets.begin(), triplets.end());

  return result;

}

VECTOR STRAND_NET_MESH::buildPerEdgeVector(const vector<VECTOR6>& perEdgeForces) const
{
  VECTOR forces(_DOFs);
  forces.setZero();

  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const VECTOR2I edge = _edgeIndices[x];
    const VECTOR6& edgeForce = perEdgeForces[x];

    unsigned int i = 3 * edge[0];
    forces[i]     += edgeForce[0];
    forces[i + 1] += edgeForce[1];
    forces[i + 2] += edgeForce[2];
    unsigned int j = 3 * edge[1];
    forces[j]     += edgeForce[3];
    forces[j + 1] += edgeForce[4];
    forces[j + 2] += edgeForce[5];
  }

  return forces;
}

SPARSE_MATRIX STRAND_NET_MESH::buildPerEdgeMatrix(const vector<MATRIX6>& perEdgeHessians) const
{
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];
    const MATRIX6& H = perEdgeHessians[i];

    // NOTE: this indexing logic has not been tested against more
    // complex geometries!
    for (int y = 0; y < 2; y++)
    {
      const int yVertex = edge[y];
      for (int x = 0; x < 2; x++)
      {
        const int xVertex = edge[x];
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

  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

// stretching 
VECTOR STRAND_NET_MESH::computeStretchingForces(const STRAND::STRETCHING& stretching) const
{
  //TIMER functionTimer(__FUNCTION__);
  vector<VECTOR6> perEdgeForces(_totalEdges);
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const VECTOR2I edge = _edgeIndices[x];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    perEdgeForces[x] = -_restEdgeLengths[x] * stretching.spatialGradient(positions, 1.0 / _restEdgeLengths[x]);
  }

  return buildPerEdgeVector(perEdgeForces);
}

SPARSE_MATRIX STRAND_NET_MESH::computeStretchingHessian(const STRAND::STRETCHING& stretching) const
{
  vector<MATRIX6> perEdgeHessians(_totalEdges);
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialHessian(positions, 1.0 / _restEdgeLengths[i]);
    perEdgeHessians[i] = -1.0 * H;
  }
  
  return buildPerEdgeMatrix(perEdgeHessians);
}

SPARSE_MATRIX STRAND_NET_MESH::computeStretchingClampedHessian(const STRAND::STRETCHING& stretching) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX6> perEdgeHessians(_totalEdges);
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialClampedHessian(positions, 1.0 / _restEdgeLengths[i]);
    perEdgeHessians[i] = -1.0 * H;
  }

  return buildPerEdgeMatrix(perEdgeHessians);
}

// bending

VECTOR STRAND_NET_MESH::computeBendingForces() 
{
  TIMER functionTimer(__FUNCTION__);

  assert(_Bs.size() == _totalBends);
  VECTOR forces(_DOFs);
  forces.setZero();
  vector<VECTOR9> perBendForces(_totalBends);
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const MATRIX2& B = _Bs[x];
    const VECTOR3I& bendVert = _bendVertices[x];
    const REAL& len = _voronoiLengths[bendVert(1)];
    BENDING_ENERGY bendMaterial(B(0, 0));
    const VECTOR3 e0 = _vertices[bendVert(1)] - _vertices[bendVert(0)];
    const VECTOR3 e1 = _vertices[bendVert(2)] - _vertices[bendVert(1)];
    MATRIX3x2 deformation; deformation << e0, e1;
    MATRIX3x2 Pk1 = bendMaterial.PK1(deformation, _restThetas[x]);
    VECTOR6 flattenPk1; flattenPk1 << Pk1.col(0), Pk1.col(1);
    const VECTOR9 force = - len * _pFpx.transpose() * flattenPk1;
    perBendForces[x] = force;
  }
  return buildPerBendVector(perBendForces);
}


SPARSE_MATRIX STRAND_NET_MESH::computeBendingHessian() 
{
  TIMER functionTimer(__FUNCTION__);
  cout<<"--------not filtered bending hessian.----------\n";
  vector<MATRIX9> perBendHessians(_totalBends);
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const MATRIX2& B = _Bs[x];
    const VECTOR3I& bendVert = _bendVertices[x];
    const REAL& len = _voronoiLengths[bendVert(1)];
    BENDING_ENERGY bendMaterial(B(0, 0));
    const VECTOR3 e0 = _vertices[bendVert(1)] - _vertices[bendVert(0)];
    const VECTOR3 e1 = _vertices[bendVert(2)] - _vertices[bendVert(1)];
    MATRIX3x2 deformation; deformation << e0, e1;
    const MATRIX9 hes = - len * _pFpx.transpose()*bendMaterial.hessian(deformation, _restThetas[x])*_pFpx;
    perBendHessians[x] = hes;
  }
  functionTimer.stop();
  return buildPerBendMatrix(perBendHessians);
}

SPARSE_MATRIX STRAND_NET_MESH::computeBendingClampedHessian()
{
  TIMER functionTimer(__FUNCTION__);
  cout<<"--------filtered bending hessian.----------\n";
  vector<MATRIX9> perBendHessians(_totalBends);
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const MATRIX2& B = _Bs[x];
    const VECTOR3I& bendVert = _bendVertices[x];
    const REAL& len = _voronoiLengths[bendVert(1)];
    BENDING_ENERGY bendMaterial(B(0, 0));
    const VECTOR3 e0 = _vertices[bendVert(1)] - _vertices[bendVert(0)];
    const VECTOR3 e1 = _vertices[bendVert(2)] - _vertices[bendVert(1)];
    MATRIX3x2 deformation; deformation << e0, e1;
    MATRIX9 hes;
    #if USING_BRUTE_FORCE_CLAMP_STRAND
    MATRIX9 hesOrig = _pFpx.transpose()*bendMaterial.hessian(deformation, _restThetas[x])*_pFpx;
    hes = - len *  clampEigenvalues(hesOrig);
    #else
    hes = - len * _pFpx.transpose()*bendMaterial.clampedHessian(deformation, _restThetas[x])*_pFpx;
    #endif
    perBendHessians[x] = hes;
  }
  functionTimer.stop();
  return buildPerBendMatrix(perBendHessians);
}

// no twisting
VECTOR STRAND_NET_MESH::computeTwistingForces() const 
{

  VECTOR forces(_DOFs);
  forces.setZero();
  return forces;
}

SPARSE_MATRIX STRAND_NET_MESH::computeTwistingHessian() const
{
  SPARSE_MATRIX hessian(_DOFs, _DOFs);
  hessian.setZero();
  return hessian;
}

SPARSE_MATRIX STRAND_NET_MESH::computeTwistingClampedHessian() const
{
  SPARSE_MATRIX hessian(_DOFs, _DOFs);
  hessian.setZero();
  return hessian;
}

} // HOBAK
