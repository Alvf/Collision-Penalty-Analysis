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
#include "STRAND_MESH.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "util/MATRIX_UTIL.h"
#include "util/TIMER.h"
#include <iostream>
#include <fstream>
#include "LINE_INTERSECT.h"
#include <float.h>

#define USING_TAN_BENDING 1
#define USING_TAN_TWISTING 1
#define USING_DWA_TWISTING 0
#define USING_GAUSS_NEWTON_TWISTING 0

Eigen::IOFormat octave(-1, 0, ", ", ";\n", "", "", "[", "]");
namespace HOBAK {

using namespace std;

// compute angle between u and v, but set the sign according to its
// projection onto n
inline REAL signedAngle(const VECTOR3& u, const VECTOR3& v, const VECTOR3& n)
{
  VECTOR3 w = u.cross(v);
  REAL angle = atan2(w.norm(), u.dot(v));
  if (n.dot(w) < 0) return -angle;
  return angle;
}

VECTOR3 findOrthogonal(const VECTOR3& u)
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

// transport vector u from t1 to t2
inline VECTOR3 parallel_transport(const VECTOR3& u,
                                  const VECTOR3& t1,
                                  const VECTOR3& t2)
{
  VECTOR3 b = t1.cross(t2);
  if (b.norm() == 0) return u;
  b.normalize();
  b = (b - (b.dot(t1) * t1)).normalized();
  b = (b - (b.dot(t2) * t2)).normalized();
  VECTOR3 n1 = t1.cross(b);
  VECTOR3 n2 = t2.cross(b);
  return u.dot(t1) * t2 + u.dot(n1) * n2 + u.dot(b) * b;
}

// rotate vector v along axis z by angle theta
// this is the Rodrigues formula: 
// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
inline VECTOR3 rotateAxisAngle(const VECTOR3& v, const VECTOR3 z, const REAL& theta)
{
  assert(fabs(z.norm() - 1.0) < 1e-10);

  if (theta == 0) return v;

  REAL c = cos(theta);
  REAL s = sin(theta);
  return c * v + s * z.cross(v) + z.dot(v) * (1.0 - c) * z;
}

///////////////////////////////////////////////////////////////////////
// assumes it's all just one big strand
///////////////////////////////////////////////////////////////////////
STRAND_MESH::STRAND_MESH(const vector<VECTOR3>& restVertices,
                         const REAL& E,        // Young's modulus
                         const REAL& G,       // shear modulus
                         const REAL& density,
                         const REAL& radiusA,
                         const REAL& radiusB) :
    _vertices(restVertices),
    _restVertices(restVertices),
    _E(E), _G(G), _density(density), _radiusA(radiusA), _radiusB(radiusB)
{
  vector<int> strand;
  for (unsigned int x = 0; x < _restVertices.size(); x++)
    strand.push_back(x);
  _strandIndices.push_back(strand);
  _totalStrands = 1;

  initialize();
  updateProperties();
}

///////////////////////////////////////////////////////////////////////
// assumes it's all just one big strand
///////////////////////////////////////////////////////////////////////
STRAND_MESH::STRAND_MESH(const vector<VECTOR3>& restVertices,
                         const vector<VECTOR3>& vertices,
                         const REAL& E,        // Young's modulus
                         const REAL& G,       // shear modulus
                         const REAL& density,
                         const REAL& radiusA,
                         const REAL& radiusB) :
    _vertices(restVertices),
    _restVertices(restVertices),
    _E(E), _G(G), _density(density), _radiusA(radiusA), _radiusB(radiusB)
{
  vector<int> strand;
  for (unsigned int x = 0; x < _restVertices.size(); x++)
    strand.push_back(x);
  _strandIndices.push_back(strand);
  _totalStrands = 1;

  initialize();

  _vertices = vertices;
  updateProperties();
  //printState();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
STRAND_MESH::STRAND_MESH(const vector<VECTOR3>& restVertices,
                         const vector<vector<int> >& strandIndices,
                         const REAL& E,        // Young's modulus
                         const REAL& G,       // shear modulus
                         const REAL& density,
                         const REAL& radiusA,
                         const REAL& radiusB) :
    _vertices(restVertices),
    _restVertices(restVertices),
    _strandIndices(strandIndices),
    _E(E), _G(G), _density(density), _radiusA(radiusA), _radiusB(radiusB)
{
  _totalStrands = _strandIndices.size();

  initialize();
  updateProperties();
  
  cout << " total strands: " << _totalStrands << endl;
  cout << " DOFs:          " << DOFs() << endl;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
STRAND_MESH::STRAND_MESH(const vector<VECTOR3>& restVertices,
                         const vector<VECTOR3>& vertices,
                         const vector<vector<int> >& strandIndices,
                         const REAL& E,        // Young's modulus
                         const REAL& G,       // shear modulus
                         const REAL& density,
                         const REAL& radiusA,
                         const REAL& radiusB) :
    _vertices(restVertices),
    _restVertices(restVertices),
    _strandIndices(strandIndices),
    _E(E), _G(G), _density(density), _radiusA(radiusA), _radiusB(radiusB)
{
  _totalStrands = _strandIndices.size();

  initialize();
  _vertices = vertices;
  updateProperties();
  
  cout << " total strands: " << _totalStrands << endl;
  cout << " DOFs:          " << DOFs() << endl;
}

///////////////////////////////////////////////////////////////////////
// generic initialization across multiple constructors
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::initialize()
{
  /*
  // initialize the edges, assuming everything is one strand
  // TODO: this needs to be made more generic
  for (unsigned int x = 1; x < _restVertices.size(); x++)
    _edgeIndices.push_back(VECTOR2I(x - 1, x));

  // initialize the bend vertices, assuming everything is one strand
  // TODO: this needs to be made more generic
  for (unsigned int x = 0; x < _restVertices.size() - 2; x++)
  {
    const VECTOR3I vertexIndices(x, x+1, x+2);
    _bendVertices.push_back(vertexIndices);
  }

  // initialize the bend edges, assuming everything is one strand
  // TODO: this needs to be made more generic
  for (unsigned int x = 0; x < _bendVertices.size(); x++)
  {
    const VECTOR2I edgeIndices(x,x+1);
    _bendEdges.push_back(edgeIndices);
  }
  */

  for (unsigned int y = 0; y < _totalStrands; y++)
  {
    const vector<int>& strand = _strandIndices[y];

    // need to store this up front, so we can index correctly
    // when making the _bendEdges
    const int edgeStart = _edgeIndices.size();

    // create the edges
    for (unsigned int x = 1; x < strand.size(); x++)
      _edgeIndices.push_back(VECTOR2I(strand[x - 1], strand[x]));

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

  /*
  for (unsigned int y = 0; y < _totalStrands; y++)
  {
    const vector<int>& strand = _strandIndices[y];
    cout << " strand " << y << ": ";
    for (unsigned int x = 0; x < strand.size(); x++)
      cout << strand[x] << " ";
    cout << endl;

    for (unsigned int x = 0; x < _edgeIndices.size(); x++)
      cout << " edge " << x << ": " << _edgeIndices[x][0] << ", " << _edgeIndices[x][1] << endl;
  }
  */

  _totalEdges = _edgeIndices.size();
  _totalBends = _bendVertices.size();
  _totalVertices = _restVertices.size();
  _DOFs = 3 * _totalVertices + _totalEdges;

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

  computeEdges();
  computeTangents();
  computeCurvatureBinormals();
  computeEdgeLengths();
  _restEdgeLengths = _edgeLengths;
  computeVoronoiLengths();

  // do the frame walk
  _director1s[0] = findOrthogonal(_tangents[0]);
  computeSpaceParallel();

  computeMaterialDirectors();
  computeVertexMasses();

  computeKappas();
  computeTwists();
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
  _bendingForceFilterEnabled = true;
  //_bendingFilterThreshold = 100.0;
  _bendingFilterThreshold = 0.0;
  //_bendingFilterThreshold = 1.0;
  
  // are we putting the edge information at the end?
  //_edgeEnd = true;
  _edgeEnd = false;

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
  buildGlobalIndices();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
STRAND_MESH::~STRAND_MESH()
{
  //if (_edgeEdgeEnergy)
  //  delete _edgeEdgeEnergy;
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setDisplacement(const VECTOR& displacements)
{
  if (_edgeEnd)
  {
    setDisplacementEdgeEnd(displacements);
    return;
  }

  setDisplacementInterleaved(displacements);
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setDisplacementInterleaved(const VECTOR& displacements)
{
  // back up old positions, for visualization
  _verticesOld = _vertices;
  _thetasOld = _thetas;

  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    int index = _globalVertexIndices[x];
    for (unsigned int y = 0; y < 3; y++, index++)
      _vertices[x][y] = _restVertices[x][y] + displacements[index];
  }
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const int index = _globalEdgeIndices[x];
    _thetas[x] = displacements[index];
  }

  updateProperties();
}

///////////////////////////////////////////////////////////////////////
// get the vertex positions, flattened into a vector
///////////////////////////////////////////////////////////////////////
const VECTOR STRAND_MESH::getDisplacement() const
{
  if (_edgeEnd)
    return getDisplacementEdgeEnd();

  return getDisplacementInterleaved();
}

///////////////////////////////////////////////////////////////////////
// get the vertex positions, flattened into a vector
///////////////////////////////////////////////////////////////////////
const VECTOR STRAND_MESH::getDisplacementInterleaved() const
{
  VECTOR result(_DOFs);
  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    const int index = _globalVertexIndices[x];
    result[index]     = _vertices[x][0] - _restVertices[x][0];
    result[index + 1] = _vertices[x][1] - _restVertices[x][1];
    result[index + 2] = _vertices[x][2] - _restVertices[x][2];
  }
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const int index = _globalEdgeIndices[x];
    result[index] = _thetas[x];
  }

  return result;
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setDisplacementEdgeEnd(const VECTOR& displacements)
{
  assert((unsigned int)displacements.size() == _DOFs);

  // back up old positions, for visualization
  _verticesOld = _vertices;
  _thetasOld = _thetas;

  unsigned int index = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++, index++)
      _vertices[x][y] = _restVertices[x][y] + displacements[index];
  }
  const unsigned int vertexEnd = 3 * _totalVertices;
  for (unsigned int x = 0; x < _totalEdges; x++)
    _thetas[x] = displacements[vertexEnd + x];

  updateProperties();
}

///////////////////////////////////////////////////////////////////////
// get the vertex positions, flattened into a vector
///////////////////////////////////////////////////////////////////////
const VECTOR STRAND_MESH::getDisplacementEdgeEnd() const
{
  VECTOR displacements(_DOFs);

  int index = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++, index++)
      displacements[index] = _vertices[x][y] - _restVertices[x][y];
  }
  const unsigned int vertexEnd = 3 * _totalVertices;
  for (unsigned int x = 0; x < _totalEdges; x++)
    displacements[vertexEnd + x] = _thetas[x];

  return displacements;
}

/*
///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setPositions(const VECTOR& positions)
{
  assert((unsigned int)positions.size() == _vertices.size() * 3);

  int index = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    for (unsigned int y = 0; y < 3; y++, index++)
      _vertices[x][y] = positions[index];

  updateProperties();
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions and thetas directly exactly
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setState(const VECTOR& state)
{
  assert(state.size() == DOFs());
  int index = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++, index++)
      _vertices[x][y] = state[index];

    index++;
    _thetas[x] = state[index];
  }

  updateProperties();
}
*/

///////////////////////////////////////////////////////////////////////
// get a vector of the vertex positions
///////////////////////////////////////////////////////////////////////
const VECTOR STRAND_MESH::getPositions() const
{
  VECTOR positions(_vertices.size() * 3);

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    positions[3 * x] = _vertices[x][0];
    positions[3 * x + 1] = _vertices[x][1];
    positions[3 * x + 2] = _vertices[x][2];
  }
  return positions;
}

/*
///////////////////////////////////////////////////////////////////////
// set the vertex positions and twist angles directly exactly
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setPositionsAndTwistAngles(const vector<VECTOR3>& positions,
                                             const VECTOR& thetas)
{
  assert(positions.size() == _vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] = positions[x];

  assert(thetas.size() == _thetas.size());
  assert(_thetas.size() == _totalEdges);
  _thetas = thetas;

  updateProperties();
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setPositions(const vector<VECTOR3>& positions)
{
  assert(positions.size() == _vertices.size());
  for (unsigned int x = 0; x < _vertices.size(); x++)
    _vertices[x] = positions[x];

  updateProperties();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setTwistAngles(const VECTOR& thetas)
{
  assert(thetas.size() == _thetas.size());
  assert(_thetas.size() == _totalEdges);
  _thetas = thetas;

  updateProperties();
}
*/

///////////////////////////////////////////////////////////////////////
// stretching energy
///////////////////////////////////////////////////////////////////////
REAL STRAND_MESH::computeStretchingEnergy(const STRAND::STRETCHING& stretching) const
{
  cout << " Stretching elements: " << endl;
  REAL totalEnergy = 0;
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const VECTOR2I& edge = _edgeIndices[x];

    // compute deformation gradient
    const VECTOR3 diff = (_vertices[edge[1]] - _vertices[edge[0]]);
    const VECTOR3 f = diff * _dmInvs[x];
    const REAL energy = stretching.psi(f);

    cout << x << ": " << energy << " dmInv: " << _dmInvs[x] << " diff.norm() " << diff.norm() << " f.norm(): " << f.norm() << endl;
    totalEnergy += _restEdgeLengths[x] * energy;
  }
  return totalEnergy;
}

///////////////////////////////////////////////////////////////////////
// Use the material PK1 to compute the stretching force
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::computeStretchingForces(const STRAND::STRETCHING& stretching) const
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

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeStretchingHessian(const STRAND::STRETCHING& stretching) const
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

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeStretchingClampedHessian(const STRAND::STRETCHING& stretching) const
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
    perEdgeHessians[i] = -1.0 * clampEigenvalues(H);
  }

  return buildPerEdgeMatrix(perEdgeHessians);
}

///////////////////////////////////////////////////////////////////////
// Use the material Hessian to compute the force gradient
// this is all super-slow, should be optimized
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeStretchingGaussNewtonHessian(const STRAND::STRETCHING& stretching) const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX6> perEdgeHessians(_totalEdges);
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    vector<VECTOR3> positions(2);
    positions[0] = _vertices[edge[0]];
    positions[1] = _vertices[edge[1]];
    const MATRIX6 H = _restEdgeLengths[i] * stretching.spatialGaussNewtonHessian(positions, 1.0 / _restEdgeLengths[i]);
    perEdgeHessians[i] = -1.0 * clampEigenvalues(H);
  }

  return buildPerEdgeMatrix(perEdgeHessians);
}

/*
///////////////////////////////////////////////////////////////////////
// change-of-basis for bending energies
// structure is:
// [ I_{3 times 3}      -I_{3 \times 3}     Zeros
// [ Zeros              -I_{3 \times 3}     I_{3 \times 3}];
///////////////////////////////////////////////////////////////////////
static MATRIX6x9 dEdxBending()
{
  MATRIX6x9 dEdx;
  dEdx.setZero();

  dEdx.block<3,3>(0,0) =  MATRIX3::Identity();
  dEdx.block<3,3>(0,3) = -MATRIX3::Identity();
  dEdx.block<3,3>(3,3) = -MATRIX3::Identity();
  dEdx.block<3,3>(3,6) =  MATRIX3::Identity();

  return dEdx;
}
*/

/*
///////////////////////////////////////////////////////////////////////
// bending energy
///////////////////////////////////////////////////////////////////////
REAL STRAND_MESH::computeBendingEnergy(const STRAND::ISOTROPIC_BENDING& bending)
{
  REAL totalEnergy = 0;
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I& bend = _bendVertices[x];
    const VECTOR3& v0 = _vertices[bend[0]];
    const VECTOR3& v1 = _vertices[bend[1]];
    const VECTOR3& v2 = _vertices[bend[2]];

    MATRIX3x2 E;
    E.col(0) = v0 - v1;
    E.col(1) = v2 - v1;

    const REAL energy = bending.psi(E, _restBendAngles[x]);
    const REAL lengthsInv = 1.0 / (_restEdgeLengths[x - 1] + _restEdgeLengths[x]);

    totalEnergy += energy * lengthsInv;
  }
  return totalEnergy;
}
*/

static MATRIX3 outerProd(const VECTOR3& x, const VECTOR3& y)
{
  return x * y.transpose();
}

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
pair<MATRIX11, MATRIX11> STRAND_MESH::computeHessianKappa(const int bendIndex) const
{
  const VECTOR2I edgeIndices = _bendEdges[bendIndex];
  const unsigned int e0 = edgeIndices[0];
  const unsigned int e1 = edgeIndices[1];

  pair<MATRIX11, MATRIX11> result;

  MATRIX11& DDkappa1 = result.first;
  MATRIX11& DDkappa2 = result.second;
  DDkappa1.setZero();
  DDkappa2.setZero();

  REAL norm_e = _edges[e0].norm();
  REAL norm_f = _edges[e1].norm();

  REAL norm2_e = norm_e * norm_e;
  REAL norm2_f = norm_f * norm_f;

  const VECTOR3& te = _tangents[e0];
  const VECTOR3& tf = _tangents[e1];

  const VECTOR3& d1e = _material1s[e0];
  const VECTOR3& d2e = _material2s[e0];
  const VECTOR3& d1f = _material1s[e1];
  const VECTOR3& d2f = _material2s[e1];

  REAL chi = 1.0 + te.dot(tf);
  VECTOR3 tilde_t = (te + tf) / chi;
  VECTOR3 tilde_d1 = (d1e + d1f) / chi;
  VECTOR3 tilde_d2 = (d2e + d2f) / chi;

  const VECTOR2& kappa = _kappas[bendIndex];
  REAL kappa1 = kappa(0);
  REAL kappa2 = kappa(1);

  const VECTOR3& kb = _kbs[bendIndex];

  /*
  std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  std::cout << " BEND " << bendIndex << std::endl;
  std::cout << " norm_e = " << norm_e << std::endl;
  std::cout << " norm_f = " << norm_e << std::endl;
  std::cout << " te = " << te.transpose().format(octave) << std::endl;
  std::cout << " tf = " << tf.transpose().format(octave) << std::endl;
  std::cout << " d1e = " << d1e.transpose().format(octave) << std::endl;
  std::cout << " d2e = " << d2e.transpose().format(octave) << std::endl;
  std::cout << " d1f = " << d1f.transpose().format(octave) << std::endl;
  std::cout << " d2f = " << d2f.transpose().format(octave) << std::endl;
  std::cout << " kappa = " << kappa.transpose().format(octave) << std::endl;
  std::cout << " kb = " << kb.transpose().format(octave) << std::endl;
  */

  MATRIX3 tt_o_tt = outerProd(tilde_t, tilde_t);
  MATRIX3 tf_c_d2t_o_tt = outerProd(tf.cross(tilde_d2), tilde_t);
  MATRIX3 tt_o_tf_c_d2t = tf_c_d2t_o_tt.transpose();
  MATRIX3 kb_o_d2e = outerProd(kb, d2e);
  MATRIX3 d2e_o_kb = kb_o_d2e.transpose();

  MATRIX3 Id = MATRIX3::Identity();

  MATRIX3 D2kappa1De2
    = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t)
    - kappa1 / (chi * norm2_e) * (Id - outerProd(te, te))
    + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

  MATRIX3 te_c_d2t_o_tt = outerProd(te.cross(tilde_d2), tilde_t);
  MATRIX3 tt_o_te_c_d2t = te_c_d2t_o_tt.transpose();
  MATRIX3 kb_o_d2f = outerProd(kb, d2f);
  MATRIX3 d2f_o_kb = kb_o_d2f.transpose();

  MATRIX3 D2kappa1Df2
    = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t)
    - kappa1 / (chi * norm2_f) * (Id - outerProd(tf, tf))
    + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

  MATRIX3 D2kappa1DeDf
    = -kappa1/(chi * norm_e * norm_f) * (Id + outerProd(te, tf))
    + 1.0 / (norm_e*norm_f) * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt + tt_o_te_c_d2t - crossProduct(tilde_d2));
  MATRIX3 D2kappa1DfDe = D2kappa1DeDf.transpose();

  MATRIX3 tf_c_d1t_o_tt = outerProd(tf.cross(tilde_d1), tilde_t);
  MATRIX3 tt_o_tf_c_d1t = tf_c_d1t_o_tt.transpose();
  MATRIX3 kb_o_d1e = outerProd(kb, d1e);
  MATRIX3 d1e_o_kb = kb_o_d1e.transpose();

  MATRIX3 D2kappa2De2
    = 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt + tt_o_tf_c_d1t)
    - kappa2 / (chi * norm2_e) * (Id - outerProd(te, te))
    - 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);

  MATRIX3 te_c_d1t_o_tt = outerProd(te.cross(tilde_d1), tilde_t);
  MATRIX3 tt_o_te_c_d1t = te_c_d1t_o_tt.transpose();
  MATRIX3 kb_o_d1f = outerProd(kb, d1f);
  MATRIX3 d1f_o_kb =  kb_o_d1f.transpose();

  MATRIX3 D2kappa2Df2
    = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - te_c_d1t_o_tt - tt_o_te_c_d1t)
    - kappa2 / (chi * norm2_f) * (Id - outerProd(tf, tf))
    - 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb);

  MATRIX3 D2kappa2DeDf
    = -kappa2/(chi * norm_e * norm_f) * (Id + outerProd(te, tf))
    + 1.0 / (norm_e*norm_f) * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t + crossProduct(tilde_d1));
  MATRIX3 D2kappa2DfDe = D2kappa2DeDf.transpose();

  REAL D2kappa1Dthetae2 = -0.5 * kb.dot(d2e);
  REAL D2kappa1Dthetaf2 = -0.5 * kb.dot(d2f);
  REAL D2kappa2Dthetae2 =  0.5 * kb.dot(d1e);
  REAL D2kappa2Dthetaf2 =  0.5 * kb.dot(d1f);

  VECTOR3 D2kappa1DeDthetae
    = 1.0 / norm_e * (0.5 * kb.dot(d1e) * tilde_t - 1.0 / chi * tf.cross(d1e));

  VECTOR3 D2kappa1DeDthetaf
    = 1.0 / norm_e * (0.5 * kb.dot(d1f) * tilde_t - 1.0 / chi * tf.cross(d1f));

  VECTOR3 D2kappa1DfDthetae
    = 1.0 / norm_f * (0.5 * kb.dot(d1e) * tilde_t + 1.0 / chi * te.cross(d1e));

  VECTOR3 D2kappa1DfDthetaf
    = 1.0 / norm_f * (0.5 * kb.dot(d1f) * tilde_t + 1.0 / chi * te.cross(d1f));

  VECTOR3 D2kappa2DeDthetae
    = 1.0 / norm_e * (0.5 * kb.dot(d2e) * tilde_t - 1.0 / chi * tf.cross(d2e));

  VECTOR3 D2kappa2DeDthetaf
    = 1.0 / norm_e * (0.5 * kb.dot(d2f) * tilde_t - 1.0 / chi * tf.cross(d2f));

  VECTOR3 D2kappa2DfDthetae
    = 1.0 / norm_f * (0.5 * kb.dot(d2e) * tilde_t + 1.0 / chi * te.cross(d2e));

  VECTOR3 D2kappa2DfDthetaf
    = 1.0 / norm_f * (0.5 * kb.dot(d2f) * tilde_t + 1.0 / chi * te.cross(d2f));

  DDkappa1.block<3,3>(0,0) =   D2kappa1De2;
  DDkappa1.block<3,3>(0,4) = - D2kappa1De2 + D2kappa1DeDf;
  DDkappa1.block<3,3>(0,8) =               - D2kappa1DeDf;
  DDkappa1.block<3,3>(4,0) = - D2kappa1De2                + D2kappa1DfDe;
  DDkappa1.block<3,3>(4,4) =   D2kappa1De2 - D2kappa1DeDf - D2kappa1DfDe + D2kappa1Df2;
  DDkappa1.block<3,3>(4,8) =                 D2kappa1DeDf                - D2kappa1Df2;
  DDkappa1.block<3,3>(8,0) =                              - D2kappa1DfDe;
  DDkappa1.block<3,3>(8,4) =                                D2kappa1DfDe - D2kappa1Df2;
  DDkappa1.block<3,3>(8,8) =                                               D2kappa1Df2;

  DDkappa1(3,3) = D2kappa1Dthetae2;
  DDkappa1(7,7) = D2kappa1Dthetaf2;

  DDkappa1.block<3,1>(0,3) = - D2kappa1DeDthetae;
  DDkappa1.block<3,1>(4,3) =   D2kappa1DeDthetae - D2kappa1DfDthetae;
  DDkappa1.block<3,1>(8,3) =                       D2kappa1DfDthetae;
  DDkappa1.block<1,3>(3,0) = DDkappa1.block<3,1>(0,3).transpose();
  DDkappa1.block<1,3>(3,4) = DDkappa1.block<3,1>(4,3).transpose();
  DDkappa1.block<1,3>(3,8) = DDkappa1.block<3,1>(8,3).transpose();

  DDkappa1.block<3,1>(0,7) = - D2kappa1DeDthetaf;
  DDkappa1.block<3,1>(4,7) =   D2kappa1DeDthetaf - D2kappa1DfDthetaf;
  DDkappa1.block<3,1>(8,7) =                       D2kappa1DfDthetaf;
  DDkappa1.block<1,3>(7,0) = DDkappa1.block<3,1>(0,7).transpose();
  DDkappa1.block<1,3>(7,4) = DDkappa1.block<3,1>(4,7).transpose();
  DDkappa1.block<1,3>(7,8) = DDkappa1.block<3,1>(8,7).transpose();

  DDkappa2.block<3,3>(0,0) =   D2kappa2De2;
  DDkappa2.block<3,3>(0,4) = - D2kappa2De2 + D2kappa2DeDf;
  DDkappa2.block<3,3>(0,8) =               - D2kappa2DeDf;
  DDkappa2.block<3,3>(4,0) = - D2kappa2De2                + D2kappa2DfDe;
  DDkappa2.block<3,3>(4,4) =   D2kappa2De2 - D2kappa2DeDf - D2kappa2DfDe + D2kappa2Df2;
  DDkappa2.block<3,3>(4,8) =                 D2kappa2DeDf                - D2kappa2Df2;
  DDkappa2.block<3,3>(8,0) =                              - D2kappa2DfDe;
  DDkappa2.block<3,3>(8,4) =                                D2kappa2DfDe - D2kappa2Df2;
  DDkappa2.block<3,3>(8,8) =                                               D2kappa2Df2;

  DDkappa2(3,3) = D2kappa2Dthetae2;
  DDkappa2(7,7) = D2kappa2Dthetaf2;

  DDkappa2.block<3,1>(0,3) = - D2kappa2DeDthetae;
  DDkappa2.block<3,1>(4,3) =   D2kappa2DeDthetae - D2kappa2DfDthetae;
  DDkappa2.block<3,1>(8,3) =                       D2kappa2DfDthetae;
  DDkappa2.block<1,3>(3,0) = DDkappa2.block<3,1>(0,3).transpose();
  DDkappa2.block<1,3>(3,4) = DDkappa2.block<3,1>(4,3).transpose();
  DDkappa2.block<1,3>(3,8) = DDkappa2.block<3,1>(8,3).transpose();

  DDkappa2.block<3,1>(0,7) = - D2kappa2DeDthetaf;
  DDkappa2.block<3,1>(4,7) =   D2kappa2DeDthetaf - D2kappa2DfDthetaf;
  DDkappa2.block<3,1>(8,7) =                       D2kappa2DfDthetaf;
  DDkappa2.block<1,3>(7,0) = DDkappa2.block<3,1>(0,7).transpose();
  DDkappa2.block<1,3>(7,4) = DDkappa2.block<3,1>(4,7).transpose();
  DDkappa2.block<1,3>(7,8) = DDkappa2.block<3,1>(8,7).transpose();

  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX11x2 STRAND_MESH::computeGradKappa(const int bendIndex) const
{
  const VECTOR2I edgeIndices = _bendEdges[bendIndex];

  const unsigned int e0 = edgeIndices[0];
  const unsigned int e1 = edgeIndices[1];

  REAL norm_e = _edges[e0].norm();
  REAL norm_f = _edges[e1].norm();

  const VECTOR3& te = _tangents[e0];
  const VECTOR3& tf = _tangents[e1];

  const VECTOR3& d1e = _material1s[e0];
  const VECTOR3& d2e = _material2s[e0];
  const VECTOR3& d1f = _material1s[e1];
  const VECTOR3& d2f = _material2s[e1];

  REAL chi = 1.0 + te.dot(tf);
  VECTOR3 tilde_t = (te + tf) / chi;
  VECTOR3 tilde_d1 = (d1e + d1f) / chi;
  VECTOR3 tilde_d2 = (d2e + d2f) / chi;

  const VECTOR2& kappa = _kappas[bendIndex];
  REAL kappa1 = kappa(0);
  REAL kappa2 = kappa(1);

  VECTOR3 Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + tf.cross(tilde_d2));
  VECTOR3 Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - te.cross(tilde_d2));

  VECTOR3 Dkappa2De = 1.0 / norm_e * (-kappa2 * tilde_t - tf.cross(tilde_d1));
  VECTOR3 Dkappa2Df = 1.0 / norm_f * (-kappa2 * tilde_t + te.cross(tilde_d1));

  MATRIX11x2 gradKappa;
  gradKappa.setZero();

  // position terms
  gradKappa.block<3,1>(0,0) = -Dkappa1De;
  gradKappa.block<3,1>(4,0) =  Dkappa1De - Dkappa1Df;
  gradKappa.block<3,1>(8,0) =              Dkappa1Df;

  gradKappa.block<3,1>(0,1) = -Dkappa2De;
  gradKappa.block<3,1>(4,1) =  Dkappa2De - Dkappa2Df;
  gradKappa.block<3,1>(8,1) =              Dkappa2Df;

  // twist terms
  const VECTOR3& kb = _kbs[bendIndex];

  gradKappa(3,0) = -0.5 * kb.dot(d1e);
  gradKappa(7,0) = -0.5 * kb.dot(d1f);
  gradKappa(3,1) = -0.5 * kb.dot(d2e);
  gradKappa(7,1) = -0.5 * kb.dot(d2f);

  return gradKappa;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL STRAND_MESH::computeBendingEnergy()
{
  assert(_Bs.size() == _totalBends);
  assert(_refVertexLengths.size() == _totalBends);
  assert(_kappas.size() == _totalBends);
  assert(_kappaBars.size() == _totalBends);

  REAL totalEnergy = 0;
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const MATRIX2& B = _Bs[x];
    const REAL& len = _refVertexLengths[x];

    const VECTOR2& kappa = _kappas[x];
    const VECTOR2& kappaBar = _kappaBars[x];

    //const REAL energy = 0.5 / len * (kappa - kappaBar).transpose() * B * (kappa - kappaBar);
    const REAL energy = 0.5 / len * (kappa - kappaBar).transpose() * B * (kappa - kappaBar);
    totalEnergy += energy;
  }

  return totalEnergy;
}

///////////////////////////////////////////////////////////////////////
// bending forces
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::computeBendingForces()
{
  return computeTanBendingForces();
}

const MATRIX6x9 dedx()
{
  MATRIX6x9 result;
  result.setZero();
  result.block<3,3>(0,0) = -MATRIX3::Identity();
  result.block<3,3>(0,3) =  MATRIX3::Identity();
  result.block<3,3>(3,3) = -MATRIX3::Identity();
  result.block<3,3>(3,6) =  MATRIX3::Identity();
  return result;
}

///////////////////////////////////////////////////////////////////////
// twisting forces
///////////////////////////////////////////////////////////////////////
MATRIX STRAND_MESH::computeTwistingHessianFiniteDiff()
{
  TIMER functionTimer(__FUNCTION__);
  MATRIX H(_DOFs, _DOFs);
  H.setZero();

  const VECTOR g0 = computeTwistingForces();
  const vector<VECTOR3> verticesOriginals = _vertices;
  const VECTOR thetasOriginals = _thetas;

  const REAL eps = 1e-8;

  //cout << " Twists: "  << _twists.transpose() << endl;
  //cout << " g0: " << g0.transpose() << endl;
  // positions first
  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    for (unsigned int y = 0; y < 3; y++)
    {
      _vertices = verticesOriginals;
      _vertices[x][y] += eps;
      computeEdges();
      computeTangents();
      computeReferenceTwist();
      //computeMaterialDirectors();
      //computeCurvatureBinormals();
      //computeKappas();

      const VECTOR g = computeTwistingForces();

      const int index = 4 * x + y;
      H.col(index) = (g - g0) / eps;
      //cout << " g: " << g.transpose() << endl;
    }
  }

  _vertices = verticesOriginals;
  _thetas = thetasOriginals;
  updateProperties();

  return H;
}

///////////////////////////////////////////////////////////////////////
// bending forces
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::computeBendingForcesFiniteDiff()
{
  TIMER functionTimer(__FUNCTION__);
  VECTOR forces(_DOFs);
  forces.setZero();

  const REAL psi0 = computeBendingEnergy();
  const vector<VECTOR3> verticesOriginals = _vertices;
  const VECTOR thetasOriginals = _thetas;

  const REAL eps = 1e-8;

  // positions first
  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    for (unsigned int y = 0; y < 3; y++)
    {
      _vertices = verticesOriginals;
      _vertices[x][y] += eps;
      computeEdges();
      computeTangents();
      computeCurvatureBinormals();
      computeKappas();

      const REAL psi = computeBendingEnergy();
      const int index = 4 * x + y;
      forces[index] = -(psi - psi0) / eps;
    }
  }
  _vertices = verticesOriginals;
  computeEdges();
  computeTangents();
  computeCurvatureBinormals();
  computeKappas();

  // twists second
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    _thetas = thetasOriginals;
    _thetas[x] += eps;
    computeReferenceTwist();
    computeMaterialDirectors();
    computeKappas();
    const REAL psi = computeBendingEnergy();
    const int index = 4 * x;
    const REAL twistForce = -(psi - psi0) / eps;
    forces[index + 3] = twistForce;
  }

  _vertices = verticesOriginals;
  _thetas = thetasOriginals;
  updateProperties();

  return forces;
}

///////////////////////////////////////////////////////////////////////
// bending forces
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::computeTanBendingForces()
{
  TIMER functionTimer(__FUNCTION__);
  assert(_Bs.size() == _totalBends);
  assert(_refVertexLengths.size() == _totalBends);
  assert(_kappas.size() == _totalBends);
  assert(_kappaBars.size() == _totalBends);

  vector<VECTOR11> perBendForces(_totalBends);

  if (_bendingForceFilterEnabled)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " BENDING FORCES FILTERED " << endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  }
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const MATRIX2& B = _Bs[x];
    const REAL& len = _refVertexLengths[x];

    const VECTOR2& kappa = _kappas[x];
    const VECTOR2& kappaBar = _kappaBars[x];

    const VECTOR11 force = -1.0/len * computeGradKappa(x) * B * (kappa - kappaBar);

    perBendForces[x] = force;

    if (_bendingForceFilterEnabled)
    {
      if (kappa.norm() > _bendingFilterThreshold)
        perBendForces[x] *= 0.0;
    }
  }

  return buildPerBendVector(perBendForces);
}

///////////////////////////////////////////////////////////////////////
// bending force gradients
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeBendingHessian()
{
  return computeTanBendingHessian();
}

///////////////////////////////////////////////////////////////////////
// bending force gradients
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTanBendingHessian()
{
  TIMER functionTimer(__FUNCTION__);
  if (_bendingForceFilterEnabled)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " BENDING FORCES FILTERED " << endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  }
  vector<MATRIX11> perBendHessians(_totalBends);
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
// bending force gradients
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeBendingClampedHessian()
{
  return computeTanBendingClampedHessian();
}

///////////////////////////////////////////////////////////////////////
// bending force gradients
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeBendingGaussNewtonHessian()
{
  return computeTanBendingGaussNewtonHessian();
}

///////////////////////////////////////////////////////////////////////
// bending force gradients
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTanBendingGaussNewtonHessian()
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX11> perBendHessians(_totalBends);

  if (_bendingForceFilterEnabled)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " BENDING FORCES FILTERED " << endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  }
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const MATRIX2& B = _Bs[x];
    const REAL len = _refVertexLengths[x];

    const VECTOR2& kappa = _kappas[x];
    const VECTOR2& kappaBar = _kappaBars[x];
    const MATRIX11x2& gradKappa = computeGradKappa(x);

    MATRIX11 localJ = 1.0 / len * gradKappa * B * gradKappa.transpose();
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
// bending force gradients
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTanBendingClampedHessian()
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX11> perBendHessians(_totalBends);

  if (_bendingForceFilterEnabled)
  {
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    std::cout << " BENDING FORCES FILTERED " << endl;
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
  }
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
VECTOR11 STRAND_MESH::computeTanGradTwist(const int bendIndex) const
{
  VECTOR11 Dtwist;
  Dtwist.setZero();

  const VECTOR2I edgeIndices = _bendEdges[bendIndex];
  const unsigned int e0 = edgeIndices[0];
  const unsigned int e1 = edgeIndices[1];

  // in case sin is being called, explicitly compute the
  // tan version here
  //VECTOR3 kb = _kbs[bendIndex];
  const VECTOR3& t0 = _tangents[e0];
  const VECTOR3& t1 = _tangents[e1];
  VECTOR3 kb = 2.0 * t0.cross(t1) / (1.0 + t0.dot(t1));

  /*
  if (_bendingForceFilterEnabled)
    if (kb.norm() > _bendingFilterThreshold)
      kb *= 0.0;
      */

  Dtwist.segment<3>(0) = -0.5 / (_edgeLengths[e0]) * kb;
  Dtwist.segment<3>(8) =  0.5 / (_edgeLengths[e1]) * kb;
  Dtwist.segment<3>(4) = -(Dtwist.segment<3>(0) + Dtwist.segment<3>(8));
  Dtwist(3) = -1;
  Dtwist(7) = 1;

  return Dtwist;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR11 STRAND_MESH::computeDWAGradTwist(const int bendIndex) const
{
  VECTOR11 Dtwist;
  Dtwist.setZero();

  const VECTOR2I edgeIndices = _bendEdges[bendIndex];
  const unsigned int i0 = edgeIndices[0];
  const unsigned int i1 = edgeIndices[1];
  const VECTOR3& e0 = _edges[i0];
  const VECTOR3& e1 = _edges[i1];
  const VECTOR3& t0 = _tangents[i0];
  const VECTOR3& t1 = _tangents[i1];
  VECTOR3 kb = (t0.cross(t1)).normalized();

  REAL angle = acos(t0.dot(t1));

  REAL coeff = tan(angle) * 0.5;  // stable, approaches the original
  //REAL coeff = angle; // also stable, doesn't resist very much
  //REAL coeff = sin(angle); // also stable, but too lively? 

  Dtwist.segment<3>(0) = -0.5 / (_edgeLengths[i0]) * coeff * kb;
  Dtwist.segment<3>(8) =  0.5 / (_edgeLengths[i1]) * coeff * kb;
  Dtwist.segment<3>(4) = -(Dtwist.segment<3>(0) + Dtwist.segment<3>(8));
  Dtwist(3) = -1;
  Dtwist(7) = 1;

  return Dtwist;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::computeTwistingForces() const
{
#if USING_DWA_TWISTING
  return computeDWATwistingForces();
#else
  return computeTanTwistingForces();
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::computeTanTwistingForces() const
{
  vector<VECTOR11> perBendForces(_totalBends);
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    REAL value = kt / len
      * (twist - undeformedTwist);

    VECTOR11 force = -value * computeTanGradTwist(x);
    /*
    VECTOR11 force = computeTanGradTwist(x);
    cout << " direction: " << force.transpose() << endl;
    cout << " magnitude: " << value << endl;
    cout <<"  twist: " << twist << " undeformed: " << undeformedTwist << endl;
    force = -value * force;
    */
    perBendForces[x] = force;
  }

  return buildPerBendVector(perBendForces);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::computeDWATwistingForces() const
{
  TIMER functionTimer(__FUNCTION__);
  vector<VECTOR11> perBendForces(_totalBends);
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    REAL value = kt / len
      * (twist - undeformedTwist);

    VECTOR11 force = -value * computeDWAGradTwist(x);
    perBendForces[x] = force;
  }

  return buildPerBendVector(perBendForces);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
MATRIX11 STRAND_MESH::computeTanHessianTwist(const int bendIndex) const
{
  MATRIX11 DDtwist;
  DDtwist.setZero();

  const VECTOR2I edgeIndices = _bendEdges[bendIndex];
  const unsigned int e0 = edgeIndices[0];
  const unsigned int e1 = edgeIndices[1];

  const VECTOR3& te = _tangents[e0];
  const VECTOR3& tf = _tangents[e1];
  const REAL norm_e = _edgeLengths[e0];
  const REAL norm_f = _edgeLengths[e1];
  VECTOR3 kb = 2.0 * te.cross(tf) / (1.0 + te.dot(tf));

  /*
  if (_bendingForceFilterEnabled)
    if (kb.norm() > _bendingFilterThreshold)
      kb *= 0.0;
      */

  REAL chi = 1 + te.dot(tf);
  VECTOR3 tilde_t = 1.0 / chi * (te + tf);

  MATRIX3 D2mDe2 = -0.25 / (norm_e * norm_e) * (outerProd(kb, te + tilde_t)
                                           + outerProd(te + tilde_t, kb));
  MATRIX3 D2mDf2 = -0.25 / (norm_f * norm_f) * (outerProd(kb, tf + tilde_t)
                                           + outerProd(tf + tilde_t, kb));
  MATRIX3 D2mDeDf = 0.5 / (norm_e * norm_f) *
                    (2.0 / chi * crossProduct(te) - outerProd(kb, tilde_t));
  MATRIX3 D2mDfDe = D2mDeDf.transpose();

  DDtwist.block<3,3>(0,0) =   D2mDe2;
  DDtwist.block<3,3>(0,4) = - D2mDe2 + D2mDeDf;
  DDtwist.block<3,3>(0,8) =          - D2mDeDf;
  DDtwist.block<3,3>(4,0) = - D2mDe2           + D2mDfDe;
  DDtwist.block<3,3>(4,4) =   D2mDe2 - D2mDeDf - D2mDfDe + D2mDf2;
  DDtwist.block<3,3>(4,8) =            D2mDeDf           - D2mDf2;
  DDtwist.block<3,3>(8,0) =                    - D2mDfDe;
  DDtwist.block<3,3>(8,4) =                      D2mDfDe - D2mDf2;
  DDtwist.block<3,3>(8,8) =                                D2mDf2;

  return DDtwist;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTwistingHessian() const
{
  TIMER functionTimer(__FUNCTION__);
#if USING_DWA_TWISTING
  return computeDWATwistingHessian();
#else
  return computeTanTwistingHessian();
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeDWATwistingClampedHessian() const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX11> perBendHessians(_totalBends);

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    const VECTOR11& gradTwist = computeDWAGradTwist(x);
    //const MATRIX11& hessTwist = computeDWAHessianTwist(x);

    MATRIX11 localJ = kt / len * (gradTwist * gradTwist.transpose());
    perBendHessians[x] = -1.0 * localJ;
  }
  
  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeDWATwistingHessian() const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX11> perBendHessians(_totalBends);

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    const VECTOR11& gradTwist = computeDWAGradTwist(x);

    MATRIX11 localJ = kt / len * (gradTwist * gradTwist.transpose());
    perBendHessians[x] = -1.0 * localJ;
  }
  
  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTanTwistingHessian() const
{
  vector<MATRIX11> perBendHessians(_totalBends);

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    const VECTOR11& gradTwist = computeTanGradTwist(x);
    const MATRIX11& hessTwist = computeTanHessianTwist(x);

#if !USING_GAUSS_NEWTON_TWISTING
    MATRIX11 localJ = kt / len * ((twist - undeformedTwist) * hessTwist
                       + gradTwist * gradTwist.transpose());
#else
    MATRIX11 localJ = kt / len * (gradTwist * gradTwist.transpose());
#endif
    perBendHessians[x] = -1.0 * localJ;

    //cout << " per bend: " << endl << perBendHessians[x] << endl;
  }
  
  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTwistingGaussNewtonHessian() const
{
  return computeTanTwistingGaussNewtonHessian();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTwistingClampedHessian() const
{
#if USING_DWA_TWISTING
  return computeDWATwistingClampedHessian();
#else
  return computeTanTwistingClampedHessian();
#endif
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTanTwistingGaussNewtonHessian() const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX11> perBendHessians(_totalBends);

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    const VECTOR11& gradTwist = computeTanGradTwist(x);
    const MATRIX11& hessTwist = computeTanHessianTwist(x);

    MATRIX11 localJ = kt / len * (gradTwist * gradTwist.transpose());
    perBendHessians[x] = -1.0 * clampEigenvalues(localJ);
  }

  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeTanTwistingClampedHessian() const
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX11> perBendHessians(_totalBends);

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const REAL len = _refVertexLengths[x];
    const REAL kt = _kts[x];
    const REAL twist = _twists[x];
    const REAL undeformedTwist = _undeformedTwists[x];

    const VECTOR11& gradTwist = computeTanGradTwist(x);
    const MATRIX11& hessTwist = computeTanHessianTwist(x);

#if !USING_GAUSS_NEWTON_TWISTING
    MATRIX11 localJ = kt / len * ((twist - undeformedTwist) * hessTwist
                       + gradTwist * gradTwist.transpose());
#else
    MATRIX11 localJ = kt / len * (gradTwist * gradTwist.transpose());
#endif
    perBendHessians[x] = -1.0 * clampEigenvalues(localJ);
  }

  return buildPerBendMatrix(perBendHessians);
}

///////////////////////////////////////////////////////////////////////
// get the area of each edge
///////////////////////////////////////////////////////////////////////
const VECTOR STRAND_MESH::voronoiAreas() const
{
  VECTOR areas(_vertices.size());
  areas.setZero();

  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    areas[x] += 0.5 * _restEdgeLengths[x];
    areas[x + 1] += 0.5 * _restEdgeLengths[x];
  }

  return areas;
}

///////////////////////////////////////////////////////////////////////
// this is here just to debug twisting energies
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::cacheOldAngles()
{
  _tangentsOld = _tangents;
  _directorOld1 = _director1s;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeEdges()
{
  // original DVT
  //for (int j = 0; j < ne(); ++j)
  //  setEdge(j, getVertex(j + 1) - getVertex(j));

  _edges.resize(_totalEdges);
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I& indices = _edgeIndices[i];
    _edges[i] = _vertices[indices[1]] - _vertices[indices[0]];

#if 1
    if (!(_edges[i].norm() > 1e-8))
    {
      cout << " vertex 0: " << _vertices[indices[0]].transpose() << endl;
      cout << " vertex 1: " << _vertices[indices[1]].transpose() << endl;

      cout << " edge indices: " << _edgeIndices[i].transpose() << endl;
      cout << " edge: " << _edges[i].transpose() << endl;
      cout << " edge norm: " << _edges[i].norm() << endl;
    }
#endif
    // edges are not degenerate, right?
    assert(_edges[i].norm() > 1e-8);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeTangents()
{
  // original DVT
  //for (int j = 0; j < ne(); ++j)
  //  setTangent(j, getEdge(j).normalized());
  assert(_tangents.size() == _edges.size());

  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    // edges are not degenerate, right?
    assert(_edges[i].norm() > 1e-8);
    _tangents[i] = _edges[i].normalized();
    
    // tangents are not degenerate, right?
    assert(_tangents[i].norm() > 1e-8);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeReferenceDirectors()
{
  assert(_director1s.size() == _totalEdges);
  assert(_director2s.size() == _totalEdges);

  // original DVT
  /*
  for (eit = edges_begin(); eit != end; ++eit) {
    edge_handle& eh = *eit;
    Vec3d t = getEdge(eh).normalized();
    Vec3d u = parallel_transport(getReferenceDirector1(eh), getTangent(eh), t);
    u = (u - u.dot(t) * t).normalized();
    setReferenceDirector1(eh, u);
    setReferenceDirector2(eh, t.cross(u));
  }
  */

  // TODO: this doesn't store an "old", so the calling order starts to matter,
  // which isn't great.
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR3 tangent = _edges[i].normalized();
    const VECTOR3 director1 = _director1s[i];
    VECTOR3 u = parallel_transport(director1, _tangents[i], tangent);
    u = (u - u.dot(tangent) * tangent).normalized();
    _director1s[i] = u;
    _director2s[i] = tangent.cross(u);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeCurvatureBinormals()
{
  computeTanBinormals();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeTanBinormals()
{
  assert(_kbs.size() == _totalBends);
  assert(_kbs.size() == _bendEdges.size());

  for (unsigned int i = 0; i < _totalBends; i++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[i];

    const VECTOR3& t0 = _tangents[edgeIndices[0]];
    const VECTOR3& t1 = _tangents[edgeIndices[1]];
    _kbs[i] = 2.0 * t0.cross(t1) / (1.0 + t0.dot(t1));
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeEdgeLengths()
{
  assert(_edgeLengths.size() == _totalEdges);

  for (unsigned int i = 0; i < _totalEdges; i++)
    _edgeLengths[i] = _edges[i].norm();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeVoronoiLengths()
{
  //setVoronoiLength(0, 0.5 * getEdgeLength(0));
  //for (int i = 1; i < nv() - 1; ++i) {
  //  setVoronoiLength(i, 0.5 * (getEdgeLength(i - 1) + getEdgeLength(i)));
  //}
  //setVoronoiLength(nv() - 1, 0.5 * getEdgeLength(ne() - 1));

  assert(_voronoiLengths.size() == _totalVertices);
  _voronoiLengths.setZero();

  // we're not tracking the endpoints, so do it this way
  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edgeIndices = _edgeIndices[i];
    const REAL& length = _edgeLengths[i];
    _voronoiLengths[edgeIndices[0]] += 0.5 * length;
    _voronoiLengths[edgeIndices[1]] += 0.5 * length;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeMaterialDirectors()
{
  assert(_material1s.size() == _totalEdges);
  assert(_material2s.size() == _totalEdges);

  for (unsigned int j = 0; j < _totalEdges; j++)
  {
    REAL c = cos(_thetas[j]);
    REAL s = sin(_thetas[j]);
    const VECTOR3& u = _director1s[j];
    const VECTOR3& v = _director2s[j];
    _material1s[j] = c * u + s * v;
    _material2s[j] = -s * u + c * v;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeReferenceTwist()
{
  //for (int i = 1; i < nv() - 1; ++i) {
  //  const Vec3d& u0 = getReferenceDirector1(i - 1);
  //  const Vec3d& u1 = getReferenceDirector1(i);
  //  const Vec3d& tangent = getTangent(i);
  //  Scalar& referenceTwist = property(m_referenceTwist)[i];

  //  // transport reference frame to next edge
  //  Vec3d ut = parallel_transport(u0, getTangent(i - 1), tangent);

  //  // rotate by current value of reference twist
  //  rotateAxisAngle(ut, tangent, referenceTwist);

  //  // compute increment to reference twist to align reference frames
  //  referenceTwist += signedAngle(ut, u1, tangent);
  //}

  assert(_referenceTwists.size() == _totalBends);

  for (unsigned int i = 0; i < _totalBends; i++)
  {
    const VECTOR2I edgeIndices = _bendEdges[i];
    const VECTOR3& u0 = _director1s[edgeIndices[0]];
    const VECTOR3& u1 = _director1s[edgeIndices[1]];
    const VECTOR3& t0 = _tangents[edgeIndices[0]];
    const VECTOR3& t1 = _tangents[edgeIndices[1]];
    const VECTOR3 uParallel = parallel_transport(u0, t0, t1);

    REAL& referenceTwist = _referenceTwists[i];
    if (t1.norm() < 1e-8)
    {
      cout << " uParallel:" << uParallel.transpose() << endl;
      cout << " t1:       " << t1.transpose() << endl;
      cout << " u0:       " << u0.transpose() << endl;
      cout << " edgeIndices:  " << edgeIndices.transpose() << endl;
      cout << " _tangents0:   " << _tangents[edgeIndices[0]].transpose() << endl;
      cout << " _tangents1:   " << _tangents[edgeIndices[1]].transpose() << endl;
    }

    // rotate by current value of reference twist
    const VECTOR3 uTransported = rotateAxisAngle(uParallel, t1, referenceTwist);

    // compute increment to reference twist to align reference frames
    referenceTwist += signedAngle(uTransported, u1, t1);
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeVertexMasses()
{
  assert(_vertexMasses.size() == _totalVertices);

  //for (unsigned int i = 0; i < _totalVertices; i++)
  //{
  //  Scalar mass = 0;
  //  if (i > 0) {
  //    mass +=
  //      computeMass(density(), radiusA(i - 1), radiusB(i - 1),
  //                  0.5 * getEdgeLength(i - 1));
  //  }
  //  if (i < nv() - 1) {
  //    mass +=
  //      computeMass(density(), radiusA(i), radiusB(i),
  //                  0.5 * getEdgeLength(i));
  //  }
  //  //assert( mass > 0.0 );
  //  setVertexMass(i, mass);
  //}

  _vertexMasses.setZero();

  for (unsigned int i = 0; i < _totalEdges; i++)
  {
    const VECTOR2I edge = _edgeIndices[i];

    //_vertexMasses[i] += density * M_PI * radiusA * radiusB * _edgeLengths[i];
    //_vertexMasses[i] += DENSITY * M_PI * RADIUSA * RADIUSB * _edgeLengths[i];
    //_vertexMasses[i] += DENSITY * M_PI * 0.001 * 0.001 * _edgeLengths[i];

    //_vertexMasses[edge[0]] += 0.5 * DENSITY * M_PI * 0.001 * 0.001 * _edgeLengths[i];
    //_vertexMasses[edge[1]] += 0.5 * DENSITY * M_PI * 0.001 * 0.001 * _edgeLengths[i];
    //_vertexMasses[edge[0]] += 0.5 * DENSITY * M_PI * RADIUSA * RADIUSB * _edgeLengths[i];
    //_vertexMasses[edge[1]] += 0.5 * DENSITY * M_PI * RADIUSA * RADIUSB * _edgeLengths[i];
    _vertexMasses[edge[0]] += 0.5 * _density * M_PI * _radiusA * _radiusB * _edgeLengths[i];
    _vertexMasses[edge[1]] += 0.5 * _density * M_PI * _radiusA * _radiusB * _edgeLengths[i];
  }
}

///////////////////////////////////////////////////////////////////////
// TODO: this needs to be made more generic, a STRAND class that
// stores each, end-to-end.
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeSpaceParallel()
{
  //// transport first edge in time
  //edge_iter eit = edges_begin(), end = edges_end();
  //edge_handle eh = *eit;
  //Vec3d t0 = getEdge(eh).normalized();
  //Vec3d u = parallel_transport(getReferenceDirector1(eh), getTangent(eh), t0);
  //u = (u - u.dot(t0) * t0).normalized();
  //setReferenceDirector1(eh, u);
  //setReferenceDirector2(eh, t0.cross(u));
  //
  //// transport along centerline (Bishop frame)
  //for (++eit; eit != end; ++eit) {
  //  eh = *eit;
  //  Vec3d t1 = getEdge(eh).normalized();
  //  u = parallel_transport(u, t0, t1);
  //  u = (u - u.dot(t1) * t1).normalized();
  //  setReferenceDirector1(eh, u);
  //  setReferenceDirector2(eh, t1.cross(u));
  //  t0 = t1;
  //}

  VECTOR3 t0 = _edges[0].normalized();
  VECTOR3 u = parallel_transport(_director1s[0], _tangents[0], t0);
  u = (u - u.dot(t0) * t0).normalized();
  _director1s[0] = u;
  _director2s[0] = t0.cross(u);

  // TODO: this needs to be generic, per strand, not ordered according to array
  for (unsigned int x = 1; x < _totalEdges; x++)
  {
    VECTOR3 t1 = _edges[x].normalized();
    u = parallel_transport(u, t0, t1);
    u = (u - u.dot(t1) * t1).normalized();
    _director1s[x] = u;
    _director2s[x] = t1.cross(u);
    t0 = t1;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeTwists()
{
  assert(_twists.size() == _totalBends);
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I edgeIndices = _bendEdges[x];

    const REAL referenceTwist = _referenceTwists[x];
    const REAL theta0 = _thetas[edgeIndices[0]];
    const REAL theta1 = _thetas[edgeIndices[1]];

    _twists[x] = referenceTwist + theta1 - theta0;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeKappas()
{
  assert(_kappas.size() == _totalBends);

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I edgeIndices = _bendEdges[x];
    const VECTOR3 kb = _kbs[x];

    const VECTOR3& m1e = _material1s[edgeIndices[0]];
    const VECTOR3& m2e = _material2s[edgeIndices[0]];
    const VECTOR3& m1f = _material1s[edgeIndices[1]];
    const VECTOR3& m2f = _material2s[edgeIndices[1]];

    MATRIX3x2 M;
    M.col(0) =  0.5 * (m2e + m2f);
    M.col(1) = -0.5 * (m1e + m1f);
    
    _kappas[x] = M.transpose() * kb;
  }
}

#if 0
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeTanKappas()
{
  assert(_kappas.size() == _totalBends);

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    //const VECTOR3I vertexIndices = _bendVertices[x];
    const VECTOR2I edgeIndices = _bendEdges[x];

    //const VECTOR3 kb = _kbs[vertexIndices[1]];
    const VECTOR3 kb = _kbs[x];

    const VECTOR3& m1e = _material1s[edgeIndices[0]];
    const VECTOR3& m2e = _material2s[edgeIndices[0]];
    const VECTOR3& m1f = _material1s[edgeIndices[1]];
    const VECTOR3& m2f = _material2s[edgeIndices[1]];

    /*
    std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    cout << " kappa before: " << _kappas[x].transpose() << endl;
    cout << " m1e: " << m1e.transpose() << endl;
    cout << " m1f: " << m1e.transpose() << endl;
    cout << " m2e: " << m2e.transpose() << endl;
    cout << " m2f: " << m2e.transpose() << endl;
    cout << " kb:  " << kb.transpose() << endl;
    */
    _kappas[x] = VECTOR2(0.5 * kb.dot(m2e + m2f), -0.5 * kb.dot(m1e + m1f));
    //cout << " kappa after: " << _kappas[x].transpose() << endl;
  }
}
#endif

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeBs(const REAL& E, const REAL& radiusA, const REAL& radiusB)
{
  assert(_Bs.size() == _totalBends);
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    MATRIX2 B;
    B.setZero();

    B(0,0) = E * M_PI * (radiusA * radiusA * radiusA * radiusB) * 0.25;
    B(1,1) = E * M_PI * (radiusA * radiusB * radiusB * radiusB) * 0.25;

    _Bs[x] = B;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeRefVertexLengths()
{
  assert(_refVertexLengths.size() == _totalBends);
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];

    const REAL length = _edgeLengths[edgeIndices[0]] + _edgeLengths[edgeIndices[1]];
    _refVertexLengths[x] = length * 0.5;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::updateProperties()
{
  computeEdges();
  computeReferenceDirectors();
  computeTangents();
  computeReferenceTwist();
  computeCurvatureBinormals();
  computeEdgeLengths();
  computeVoronoiLengths();
  computeMaterialDirectors();

  // Updating so that density update is taken into account
  computeVertexMasses();

  computeKappas();
  computeTwists();
}

///////////////////////////////////////////////////////////////////////
// print out all the vertex positions and thetas
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::printState()
{
  std::cout << "==================================================" << std::endl;
  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    std::cout << " vertex " << x << ": " << _vertices[x].transpose() << std::endl;
  }
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    std::cout << " theta " << x << ": " << _thetas[x] << std::endl;
  }
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    //const Scalar referenceTwist = this->getReferenceTwist(x);
    std::cout << " reference twist " << x << ": " << _referenceTwists[x] << std::endl;
  }
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    std::cout << " curvature binormal " << x << ": " << _kbs[x].transpose() << std::endl;
  }
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    std::cout << " reference 1 director " << x << ": " << _director1s[x].transpose() << std::endl;
  }
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    std::cout << " reference 2 director " << x << ": " << _director2s[x].transpose() << std::endl;
  }
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    std::cout << " material 1 director " << x << ": " << _material1s[x].transpose() << std::endl;
  }
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    std::cout << " material 2 director " << x << ": " << _material2s[x].transpose() << std::endl;
  }
  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    std::cout << " mass " << x << ": " << _vertexMasses[x] << std::endl;
  }
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    std::cout << " edge length " << x << ": " << _edgeLengths[x] << std::endl;
  }
  std::cout << "==================================================" << std::endl;
}
/*
{
  cout << " Vertices: " << endl;
  for (unsigned int x = 0; x < _vertices.size(); x++)
    cout << _vertices[x].transpose() << endl;
  cout << " Thetas: " << endl;
  for (unsigned int x = 0; x < _thetas.size(); x++)
    cout << _thetas[x] << endl;
}
*/

///////////////////////////////////////////////////////////////////////
// edge-edge collision detection
//
// TODO: each edge only owns one vertex; the downstream one
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::computeEdgeEdgeCollisions(const REAL& collisionEps, const bool verbose)
{
  TIMER functionTimer(__FUNCTION__);

  //const REAL collisionEps = edgeEdgeEnergy.eps();

  // for visualization purposes
  _edgeEdgeCollisionsOld = _edgeEdgeCollisions;
  _edgeEdgeCoordinatesOld = _edgeEdgeCoordinates;

  _edgeEdgeCollisions.clear();
  _edgeEdgeIntersections.clear();
  _edgeEdgeCoordinates.clear();
  _edgeEdgeCollisionAreas.clear();

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
#if 0
    _collisionTree->nearbyEdges(_edgeIndices[x], _collisionEps, nearbyEdges);
#else
    for (unsigned int y = 0; y < _edgeIndices.size(); y++)
      nearbyEdges.push_back(y);
#endif

    // find the closest other edge
    for (unsigned int y = 0; y < nearbyEdges.size(); y++)
    {
      if (x == y) continue;

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

      /*
      const REAL innerDiff = (middle0 - outerPoint).norm();
      const REAL outerDiff = (middle1 - innerPoint).norm();
      cout << " middle diffs:" << outerDiff << " " << innerDiff << endl;
      */

      //if (diff.norm() > _collisionEps) continue;
      if (diff.norm() > collisionEps) continue;

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
#if 0
{
  TIMER functionTimer(__FUNCTION__);

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
    int closestEdge = -1;
    REAL closestDistance = FLT_MAX;
    VECTOR2 aClosest(-1,-1);
    VECTOR2 bClosest(-1,-1);
    const VECTOR2I& outerEdge = _edgeIndices[x];
    const VECTOR3& v0 = _vertices[outerEdge[0]];
    const VECTOR3& v1 = _vertices[outerEdge[1]];
    //const unsigned int outerFlat = outerEdge[0] + outerEdge[1] * _edgeIndices.size();

    vector<int> nearbyEdges;
#if 0
    _collisionTree->nearbyEdges(_edgeIndices[x], _collisionEps, nearbyEdges);
#else
    for (unsigned int y = 0; y < _edgeIndices.size(); y++)
      nearbyEdges.push_back(y);
#endif

    // find the closest other edge
    for (unsigned int y = 0; y < nearbyEdges.size(); y++)
    {
      if (x == y) continue;

      /*
      bool debug = false;
      //if (x == 2 && nearbyEdges[y] == 5)
      if (x == 2)
      {
        debug = true;
        std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
        cout << " HERE " << endl;
      }
      */

      // skip if index is smaller -- don't want to double count nearby edges
      // (a,b) and (b,a)
      if ((unsigned int)nearbyEdges[y] < x) continue;

      const VECTOR2I innerEdge = _edgeIndices[nearbyEdges[y]];
      /*
      if (debug)
      {
        std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
        cout << " inner: " << innerEdge.transpose() << endl;
        cout << " outer: " << outerEdge.transpose() << endl;
      }
      */

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
      //const REAL skipEps = 1e-4;
      const REAL skipEps = 0;
      if ((a[0] < skipEps) || (a[0] > 1.0 - skipEps)) continue;
      if ((a[1] < skipEps) || (a[1] > 1.0 - skipEps)) continue;
      if ((b[0] < skipEps) || (b[0] > 1.0 - skipEps)) continue;
      if ((b[1] < skipEps) || (b[1] > 1.0 - skipEps)) continue;

      // it's mid-segment, and closest, so remember it
      closestDistance = distance;
      closestEdge = nearbyEdges[y];

      aClosest = a;
      bClosest = b;
    }

    //std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
    //cout << " Closest edge: " << closestDistance << endl;

    // if nothing was close, move on
    if (closestEdge == -1) continue;

    // if it's within the positive threshold, it's in collision
    if (closestDistance < _collisionEps)
    {
      pair<int,int> collision(x, closestEdge);
      _edgeEdgeCollisions.push_back(collision);

      // this was actually set, right?
      //assert(aClosest[0] > 0.0 && aClosest[1] > 0.0);
      //assert(bClosest[0] > 0.0 && bClosest[1] > 0.0);

      pair<VECTOR2,VECTOR2> coordinate(aClosest, bClosest);
      _edgeEdgeCoordinates.push_back(coordinate);

      // get the areas too
      const VECTOR2I innerEdge = _edgeIndices[closestEdge];
      const pair<int,int> outerPair(outerEdge[0], outerEdge[1]);
      const pair<int,int> innerPair(innerEdge[0], innerEdge[1]);
      const REAL xArea = _restEdgeLengths[_edgeHash[outerPair]];
      const REAL closestArea = _restEdgeLengths[_edgeHash[innerPair]];
      _edgeEdgeCollisionAreas.push_back(xArea + closestArea);
    }
  }
  assert(_edgeEdgeCollisions.size() == _edgeEdgeCoordinates.size());

  //if (_edgeEdgeCollisions.size() > 0)
    cout << " Found " << _edgeEdgeCollisions.size() << " edge-edge collisions " << endl;

//#define VERY_VERBOSE 1
#if 1
  if (_edgeEdgeCollisions.size() > 0)
    cout << " pairs: " << endl;
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
#endif

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision forces using x-based formulation
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::computeEdgeEdgeCollisionForces(const ENERGY_12D& eeEnergy) const
{
  TIMER functionTimer(__FUNCTION__);

  vector<VECTOR12> perEdgeForces(_edgeEdgeCollisions.size());
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edgeIndices[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edgeIndices[_edgeEdgeCollisions[i].second];

    vector<VECTOR3> vs(4);
    vs[0] = _vertices[edge0[0]];
    vs[1] = _vertices[edge0[1]];
    vs[2] = _vertices[edge1[0]];
    vs[3] = _vertices[edge1[1]];
    VECTOR3 e0n = (vs[1] - vs[0]).normalized();
    VECTOR3 e1n = (vs[3] - vs[2]).normalized();
    if (abs(e0n.dot(e1n)) > 1 - 1e-2) {//near-parallel guard
      perEdgeForces[i] = VECTOR12::Zero();
      continue;  
    }
    
    vector<int> extraInfo;
    extraInfo.push_back(1);
    REAL len = eeEnergy.getLen(eeEnergy.makeState(vs));
    if(len < 0){
      extraInfo[0] = -1;
    }

    const VECTOR12 force = -_edgeEdgeCollisionAreas[i] * eeEnergy.gradient(vs,extraInfo);

    perEdgeForces[i] = force;
  }

  return buildEdgeEdgeVector(perEdgeForces);
}

///////////////////////////////////////////////////////////////////////
// compute edge-edge collision Hessians using x-based formulation
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::computeEdgeEdgeCollisionClampedHessian(const ENERGY_12D& eeEnergy) const
{
  TIMER functionTimer(__FUNCTION__);

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
    VECTOR3 e0n = (vs[1] - vs[0]).normalized();
    VECTOR3 e1n = (vs[3] - vs[2]).normalized();
    if (abs(e0n.dot(e1n)) > 1 - 1e-2) {//near-parallel guard
      perEdgeHessians[i] = MATRIX12::Zero();
      continue;  
    }
    
    vector<int> extraInfo;
    extraInfo.push_back(1);
    REAL len = eeEnergy.getLen(eeEnergy.makeState(vs));
    if(len < 0){
      extraInfo[0] = -1;
    }

    const MATRIX12 H = -_edgeEdgeCollisionAreas[i] * eeEnergy.clampedHessian(vs,extraInfo);
    perEdgeHessians[i] = H;
  }

  return buildEdgeEdgeMatrix(perEdgeHessians);
}

///////////////////////////////////////////////////////////////////////
// find the strand with the longest length
///////////////////////////////////////////////////////////////////////
const int STRAND_MESH::longestStrand() const
{
  int longestIndex = 0;
  REAL longestLength = strandLength(0);

  for (unsigned int x = 0; x < _strandIndices.size(); x++)
  {
    const REAL length = strandLength(x);

    if (length > longestLength)
    {
      longestIndex = x;
      longestLength = length;
    }
  }
  return longestIndex;
}

///////////////////////////////////////////////////////////////////////
// get the length of an individual strand
///////////////////////////////////////////////////////////////////////
const REAL STRAND_MESH::strandLength(const unsigned int index) const
{
  const vector<int>& strand = _strandIndices[index];

  REAL length = 0;
  for (unsigned int x = 0; x < strand.size() - 1; x++)
  {
    const VECTOR3& v0 = _vertices[strand[x]];
    const VECTOR3& v1 = _vertices[strand[x + 1]];

    length += (v0 - v1).norm();
  }

  return length;
}

///////////////////////////////////////////////////////////////////////
// output a single strand to an SOBJ file
///////////////////////////////////////////////////////////////////////
const bool STRAND_MESH::writeStrand(const string& filename, const int index) const
{
  FILE* file = fopen(filename.c_str(), "w");
  if (file == NULL)
  {
    cout << " Failed to open file " << filename.c_str() << endl;
    return false;
  }

  if (index >= (int)_strandIndices.size()) return false;
 
  // get the strand 
  const vector<int>& strand = _strandIndices[index];

  // write out its vertices
  for (unsigned int x = 0; x < strand.size(); x++)
  {
    const VECTOR3& v = _vertices[strand[x]];
    fprintf(file, "v %f %f %f\n", v[0], v[1], v[2]);
  }

  // the strand index is then just canonical
  fprintf(file, "s ");
  for (unsigned int x = 0; x < strand.size(); x++)
    fprintf(file, "%i ", x);
  fprintf(file, "\n\n");
  
  fclose(file);
  cout << " done. " << endl;
  return true;
}

///////////////////////////////////////////////////////////////////////
// shared matrix construction between clamped and unclamped
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildPerBendMatrix(const vector<MATRIX11>& perBendHessians) const
{
  if (_edgeEnd)
    return buildPerBendMatrixEdgeEnd(perBendHessians);
  return buildPerBendMatrixInterleaved(perBendHessians);
}
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildPerEdgeVector(const vector<VECTOR6>& perEdgeForces) const
{
  if (_edgeEnd)
    return buildPerEdgeVectorEdgeEnd(perEdgeForces);
  return buildPerEdgeVectorInterleaved(perEdgeForces);
}

///////////////////////////////////////////////////////////////////////
// distribute the edge Hessians to the global matrix
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildPerEdgeMatrix(const vector<MATRIX6>& perEdgeHessians) const
{
  if (_edgeEnd)
    return buildPerEdgeMatrixEdgeEnd(perEdgeHessians);
  return buildPerEdgeMatrixInterleaved(perEdgeHessians);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildPerBendVector(const vector<VECTOR11>& perBendForces) const
{
  if (_edgeEnd)
    return buildPerBendVectorEdgeEnd(perBendForces);
  return buildPerBendVectorInterleaved(perBendForces);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildEdgeEdgeVector(const vector<VECTOR12>& perEdgeForces) const
{
  if (_edgeEnd)
    return buildEdgeEdgeVectorEdgeEnd(perEdgeForces);
  return buildEdgeEdgeVectorInterleaved(perEdgeForces);
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildEdgeEdgeMatrix(const vector<MATRIX12>& perEdgeHessians) const
{
  if (_edgeEnd)
    return buildEdgeEdgeMatrixEdgeEnd(perEdgeHessians);
  return buildEdgeEdgeMatrixInterleaved(perEdgeHessians);
}

///////////////////////////////////////////////////////////////////////
// shared matrix construction between clamped and unclamped
// assumes that all the edge DOFs are packed at the end
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildPerBendMatrixEdgeEnd(const vector<MATRIX11>& perBendHessians) const
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

///////////////////////////////////////////////////////////////////////
// shared matrix construction between clamped and unclamped
// assumes that all the edge DOFs are interleaved after vertices
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildPerBendMatrixInterleaved(const vector<MATRIX11>& perBendHessians) const
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
      const unsigned int viGlobal = _globalVertexIndices[vi];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        const unsigned int vjGlobal = _globalVertexIndices[vj];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
           {
              const REAL entry = H(4 * i + a, 4 * j + b);
              TRIPLET triplet(viGlobal + a, vjGlobal + b,entry);
              triplets.push_back(triplet);
           }
      }
    }
  }

  // do the edge-edge (twist) DOFs second
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];
    const MATRIX11& H = perBendHessians[x];

    for (unsigned int i = 0; i < 2; i++)
    {
      const unsigned int ei = edgeIndices[i];
      const unsigned int eiGlobal = _globalEdgeIndices[ei];
      for (unsigned int j = 0; j < 2; j++)
      {
        const unsigned int ej = edgeIndices[j];
        const unsigned int ejGlobal = _globalEdgeIndices[ej];
        TRIPLET triplet(eiGlobal, ejGlobal, H(4 * i + 3, 4 * j + 3));
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
      const unsigned int eiGlobal = _globalEdgeIndices[ei];
      for (unsigned int j = 0; j < 3; j++)
      {
        const unsigned int vj = vertexIndices[j];
        const unsigned int vjGlobal = _globalVertexIndices[vj];
        for (unsigned int a = 0; a < 3; a++)
        {
          const unsigned int row = eiGlobal;
          const unsigned int col = vjGlobal + a;
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

///////////////////////////////////////////////////////////////////////
// assumes that all the edge DOFs are packed at the end
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildPerEdgeVectorEdgeEnd(const vector<VECTOR6>& perEdgeForces) const
{
  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3 + _totalEdges;
  VECTOR forces(DOFs);
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

///////////////////////////////////////////////////////////////////////
// assumes that all the edge DOFs are interleaved with vertex DOFs
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildPerEdgeVectorInterleaved(const vector<VECTOR6>& perEdgeForces) const
{
  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3 + _totalEdges;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const VECTOR2I edge = _edgeIndices[x];
    const VECTOR6& edgeForce = perEdgeForces[x];

    unsigned int i = _globalVertexIndices[edge[0]];
    forces[i]     += edgeForce[0];
    forces[i + 1] += edgeForce[1];
    forces[i + 2] += edgeForce[2];
    unsigned int j = _globalVertexIndices[edge[1]];
    forces[j]     += edgeForce[3];
    forces[j + 1] += edgeForce[4];
    forces[j + 2] += edgeForce[5];
  }

  return forces;
}

///////////////////////////////////////////////////////////////////////
// distribute the edge Hessians to the global matrix
// assumes that all the edge DOFs are packed at the end
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildPerEdgeMatrixEdgeEnd(const vector<MATRIX6>& perEdgeHessians) const
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

  int DOFs = _vertices.size() * 3 + _totalEdges;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// distribute the edge Hessians to the global matrix
// assumes that all the edge DOFs are interleaved with vertices
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildPerEdgeMatrixInterleaved(const vector<MATRIX6>& perEdgeHessians) const
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
      const int yGlobal = _globalVertexIndices[yVertex];
      for (int x = 0; x < 2; x++)
      {
        const int xVertex = edge[x];
        const int xGlobal = _globalVertexIndices[xVertex];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            const REAL entry = H(3 * x + a, 3 * y + b);
            TRIPLET triplet(xGlobal + a, yGlobal + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  int DOFs = _vertices.size() * 3 + _totalEdges;
  SPARSE_MATRIX A(DOFs, DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// assumes that all the edge DOFs are packed at the end
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildPerBendVectorEdgeEnd(const vector<VECTOR11>& perBendForces) const 
{
  VECTOR forces(_DOFs);
  forces.setZero();
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I vertexIndices = _bendVertices[x];
    const unsigned int v0 = vertexIndices[0];
    const unsigned int v1 = vertexIndices[1];
    const unsigned int v2 = vertexIndices[2];

    forces[3 * v0]     += perBendForces[x][0];
    forces[3 * v0 + 1] += perBendForces[x][1];
    forces[3 * v0 + 2] += perBendForces[x][2];

    forces[3 * v1]     += perBendForces[x][4];
    forces[3 * v1 + 1] += perBendForces[x][5];
    forces[3 * v1 + 2] += perBendForces[x][6];

    forces[3 * v2]     += perBendForces[x][8];
    forces[3 * v2 + 1] += perBendForces[x][9];
    forces[3 * v2 + 2] += perBendForces[x][10];
  }

  const int vertexEnd = 3 * _totalVertices;
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];
    const unsigned int e0 = edgeIndices[0];
    const unsigned int e1 = edgeIndices[1];

    forces[vertexEnd + e0] += perBendForces[x][3];
    forces[vertexEnd + e1] += perBendForces[x][7];
  }

  return forces;
}

///////////////////////////////////////////////////////////////////////
// assumes that all the edge DOFs are interleaved with vertex DOFs
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildPerBendVectorInterleaved(const vector<VECTOR11>& perBendForces) const 
{
  VECTOR forces(_DOFs);
  forces.setZero();
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I vertexIndices = _bendVertices[x];
    const unsigned int v0 = _globalVertexIndices[vertexIndices[0]];
    const unsigned int v1 = _globalVertexIndices[vertexIndices[1]];
    const unsigned int v2 = _globalVertexIndices[vertexIndices[2]];

    forces[v0]     += perBendForces[x][0];
    forces[v0 + 1] += perBendForces[x][1];
    forces[v0 + 2] += perBendForces[x][2];

    forces[v1]     += perBendForces[x][4];
    forces[v1 + 1] += perBendForces[x][5];
    forces[v1 + 2] += perBendForces[x][6];

    forces[v2]     += perBendForces[x][8];
    forces[v2 + 1] += perBendForces[x][9];
    forces[v2 + 2] += perBendForces[x][10];
  }

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR2I& edgeIndices = _bendEdges[x];
    const unsigned int e0 = _globalEdgeIndices[edgeIndices[0]];
    const unsigned int e1 = _globalEdgeIndices[edgeIndices[1]];

    forces[e0] += perBendForces[x][3];
    forces[e1] += perBendForces[x][7];
  }

  return forces;
}

///////////////////////////////////////////////////////////////////////
// assumes that all the edge DOFs are packed at the end
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildEdgeEdgeVectorEdgeEnd(const vector<VECTOR12>& perEdgeForces) const
{
  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  VECTOR forces(_DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edgeIndices[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edgeIndices[_edgeEdgeCollisions[i].second];
    const VECTOR12& edgeForce = perEdgeForces[i];

    vector<int> vertexIndices(4);
    vertexIndices[0] = edge0[0];
    vertexIndices[1] = edge0[1];
    vertexIndices[2] = edge1[0];
    vertexIndices[3] = edge1[1];

    for (int x = 0; x < 4; x++)
    {
      unsigned int index = 3 * vertexIndices[x];
      assert(index < _DOFs);
      forces[index]     += edgeForce[3 * x];
      forces[index + 1] += edgeForce[3 * x + 1];
      forces[index + 2] += edgeForce[3 * x + 2];
    }
  }

  return forces;
}

///////////////////////////////////////////////////////////////////////
// assumes that all the edge DOFs are interleaved with vertex DOFs
///////////////////////////////////////////////////////////////////////
VECTOR STRAND_MESH::buildEdgeEdgeVectorInterleaved(const vector<VECTOR12>& perEdgeForces) const
{
  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  VECTOR forces(_DOFs);
  forces.setZero();

  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const VECTOR2I& edge0 = _edgeIndices[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edgeIndices[_edgeEdgeCollisions[i].second];
    const VECTOR12& edgeForce = perEdgeForces[i];

    vector<int> vertexIndices(4);
    vertexIndices[0] = _globalVertexIndices[edge0[0]];
    vertexIndices[1] = _globalVertexIndices[edge0[1]];
    vertexIndices[2] = _globalVertexIndices[edge1[0]];
    vertexIndices[3] = _globalVertexIndices[edge1[1]];

    for (int x = 0; x < 4; x++)
    {
      unsigned int index = vertexIndices[x];
      assert(index < _DOFs);
      forces[index]     += edgeForce[3 * x];
      forces[index + 1] += edgeForce[3 * x + 1];
      forces[index + 2] += edgeForce[3 * x + 2];
    }
  }

  return forces;
}

///////////////////////////////////////////////////////////////////////
// assumes that all the edge DOFs are packed at the end
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildEdgeEdgeMatrixEdgeEnd(const vector<MATRIX12>& perEdgeHessians) const
{
  TIMER functionTimer(__FUNCTION__);
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const MATRIX12& H = perEdgeHessians[i];
    const VECTOR2I& edge0 = _edgeIndices[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edgeIndices[_edgeEdgeCollisions[i].second];

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

  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}

///////////////////////////////////////////////////////////////////////
// assumes that all the edge DOFs are interleaved with vertex DOFs
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX STRAND_MESH::buildEdgeEdgeMatrixInterleaved(const vector<MATRIX12>& perEdgeHessians) const
{
  TIMER functionTimer(__FUNCTION__);
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const MATRIX12& H = perEdgeHessians[i];
    const VECTOR2I& edge0 = _edgeIndices[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edgeIndices[_edgeEdgeCollisions[i].second];

    vector<int> vertexIndices(4);
    vertexIndices[0] = _globalVertexIndices[edge0[0]];
    vertexIndices[1] = _globalVertexIndices[edge0[1]];
    vertexIndices[2] = _globalVertexIndices[edge1[0]];
    vertexIndices[3] = _globalVertexIndices[edge1[1]];

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
            TRIPLET triplet(xVertex + a, yVertex + b, entry);
            triplets.push_back(triplet);
          }
      }
    }
  }

  SPARSE_MATRIX A(_DOFs, _DOFs);
  A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}
  
///////////////////////////////////////////////////////////////////////
// build out the indexing into the global vector
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::buildGlobalIndices()
{
  // see which edge each vertex owns
  _vertexOwnsEdge.resize(_totalVertices);
  for (unsigned int x = 0; x < _totalVertices; x++)
    _vertexOwnsEdge[x] = -1;

  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I vertices = _bendVertices[x];
    const VECTOR2I edges = _bendEdges[x];

    // this one actually may have been assigned in
    // the previous bend
    // assert(_vertexOwnsEdge[vertices[0]] == -1);
    
    // this hasn't been assigned before, right?
    assert(_vertexOwnsEdge[vertices[1]] == -1);

    _vertexOwnsEdge[vertices[0]] = edges[0];
    _vertexOwnsEdge[vertices[1]] = edges[1];
  }

  _globalVertexIndices.resize(_totalVertices);
  _globalEdgeIndices.resize(_totalEdges);
  _globalStrandEnds.resize(_totalStrands);

  unsigned int currentStrand = 0;
  unsigned int index = 0;
  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    _globalVertexIndices[x] = index;
    index += 3;

    // if it owns an edge, store it
    if (_vertexOwnsEdge[x] != -1)
    {
      const int whichEdge = _vertexOwnsEdge[x];
      _globalEdgeIndices[whichEdge] = index;
      index++;
    }
    // if it doesn't own an edge, it must be the end of a strand
    // so store that as well
    else
    {
      // DEBUG: why is this guard needed for DVT?
      if (currentStrand < _totalStrands)
        _globalStrandEnds[currentStrand] = index;
      currentStrand++;
    }
  }

  // this all lines up, right?
  assert(index == _DOFs);
  //assert(currentStrand == _totalStrands);
}

///////////////////////////////////////////////////////////////////////
// for scene replay visualization, allow collision pairs to be reset
///////////////////////////////////////////////////////////////////////
void STRAND_MESH::setCollisionPairs(const vector<pair<int,int>>& eeColllsions)
{
  _edgeEdgeCollisions = eeColllsions;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool STRAND_MESH::readSOBJFile(const char *filename, vector<VECTOR3> &vertices,
                               vector<vector<int>> &indices)
{
  vertices.clear();
  indices.clear();
  std::ifstream ifs(filename);
  if (!ifs.good()) {
    std::cerr << "Cannot read file at path " << filename << std::endl;
    return false;
  }

  const std::string pound("#");
  const std::string v("v");
  const std::string l("l");
  const std::string s("s");

  for (std::string line; std::getline(ifs, line);) {
    std::istringstream iss(line);
    std::string key;
    iss >> key;

    if (key == v) {
      double x, y, z;

      // Get the vertex positions.
      iss >> x >> y >> z;
      vertices.emplace_back(x, y, z);
    } else if (key == l) {
      // Ignore this key right now since we want it by strand.
    } else if (key == s) {
      // Add a new index vector
      indices.emplace_back(vector<int>{});
      std::string edge;
      while (!iss.eof()) {
        // End of the road
        if (iss.tellg() == -1) {
          break;
        }
        // Fetch the index otherwise
        iss >> edge;
        // Load into the index vector
        indices.at(indices.size() - 1).push_back(std::stoi(edge));

        //cout << " edge: " << edge << endl;
      }
    } else if (key == pound) {
      // Skip comments
      continue;
    } else {
      std::cerr << "Found invalid line in sobj " << line << std::endl;
    }
  }

  return true;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
bool STRAND_MESH::writeSOBJFile(const char *filename, const vector<VECTOR3> &vertices,
                                const vector<vector<int>> &indices)
{
  // Now, write the points and indices
  std::ofstream s(filename);
  if (!s.is_open()) {
    std::cerr << "IOError: write_strand_point_mesh() could not open "
              << filename << std::endl;
    return false;
  }

  if (vertices.empty()) {
    std::cerr << "No vertices to write" << std::endl;
    return false;
  }

  if (indices.empty()) {
    std::cerr << "No indices to write" << std::endl;
    return false;
  }

  for (int ii = 0; ii < vertices.size(); ++ii) {
    const VECTOR3 &row = vertices.at(ii);
    const auto x = row.x();
    const auto y = row.y();
    const auto z = row.z();

    s << "v " << x << " " << y << " " << z << "\n";
  }

  for (int ii = 0; ii < indices.size(); ++ii) {
    // We only want the single index since they just duplicate anyway
    auto dedupedIndices = indices.at(ii);
    dedupedIndices.erase(
        std::unique(dedupedIndices.begin(), dedupedIndices.end()),
        dedupedIndices.end());
    s << "s ";
    for (int jj = 0; jj < dedupedIndices.size() - 1; ++jj) {
      const auto& index = dedupedIndices.at(jj);
      s << index << " ";
    }
    s << dedupedIndices.at(dedupedIndices.size() - 1) << "\n";
  }

  return true;
}


}
