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
#include "TET_STRAND_MESH.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_UNIT_BENDING.h"
#include "Hyperelastic/Strand/QUADRATIC_STRETCHING.h"
#include "Hyperelastic/Volume/ANISOTROPIC_ARAP.h"
#include "Hyperelastic/Volume/ANISOTROPIC_DIRICHLET.h"
#include "Hyperelastic/Volume/TET_STRAND_TWIST.h"
#include "util/MATRIX_UTIL.h"
#include "util/TIMER.h"
#include <iostream>
#include "LINE_INTERSECT.h"
#include <float.h>
#include "Hyperelastic/Volume/ARAP.h"

namespace HOBAK {

using namespace std;

///////////////////////////////////////////////////////////////////////
// accepts a vector of individual strands
///////////////////////////////////////////////////////////////////////
TET_STRAND_MESH::TET_STRAND_MESH(const vector<VECTOR3>& restVertices,
                         const vector<vector<int> >& strandIndices,
                         const REAL& E,        // Young's modulus
                         const REAL& G,       // shear modulus
                         const REAL& density,
                         const REAL& radiusA,
                         const REAL& radiusB) :
  STRAND_MESH(restVertices, strandIndices, E, G, density, radiusA, radiusB)
{
  _DOFs = 3 * _vertices.size();
  _edgeEnd = true;
  initializeTets();
  initializeFastHessian();
 
  // rebuild global indices without an edge twist
  rebuildGlobalIndices();

  _collisionTree = new AABB_TREE(_vertices, &_edgeIndices);

  _strandSelfCollisionDisabled = false;
  buildEdgeStrandIndices();
}

///////////////////////////////////////////////////////////////////////
// accepts a vector of individual strands, and a modified initial pose
///////////////////////////////////////////////////////////////////////
TET_STRAND_MESH::TET_STRAND_MESH(const vector<VECTOR3>& restVertices,
                         const vector<VECTOR3>& vertices,
                         const vector<vector<int> >& strandIndices,
                         const REAL& E,        // Young's modulus
                         const REAL& G,       // shear modulus
                         const REAL& density,
                         const REAL& radiusA,
                         const REAL& radiusB) :
  STRAND_MESH(restVertices, vertices, strandIndices, E, G, density, radiusA, radiusB)
{
  _DOFs = 3 * _vertices.size();
  _edgeEnd = true;
  initializeTets();
  initializeFastHessian();
  
  // rebuild global indices without an edge twist
  rebuildGlobalIndices();

  _collisionTree = new AABB_TREE(_vertices, &_edgeIndices);
  
  _strandSelfCollisionDisabled = false;
  buildEdgeStrandIndices();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
TET_STRAND_MESH::~TET_STRAND_MESH()
{
  if (_collisionTree)
    delete _collisionTree;
}

///////////////////////////////////////////////////////////////////////
// find which strand each edge belongs to
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::buildEdgeStrandIndices()
{
  _perEdgeStrandIndex.resize(_totalEdges);

  for (unsigned int x = 0; x < _strandIndices.size(); x++)
  {
    for (unsigned int y = 1; y < _strandIndices[x].size(); y++)
    {
      // find the vertex index pair
      const pair<int,int> edgePair(_strandIndices[x][y - 1], _strandIndices[x][y]);
      assert(_edgeHash.find(edgePair) != _edgeHash.end());

      // find the edge index in the hash
      const int edgeIndex = _edgeHash[edgePair];
      assert(edgeIndex < _totalEdges);

      // store the strand index
      _perEdgeStrandIndex[edgeIndex] = x;
    }
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::initializeFastHessian()
{
  cout << " Initializing matrix sparsity ... " << flush;
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  triplets.reserve(_tets.size() * 144);
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    const MATRIX12& H = MATRIX12::Zero();
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

  // bake out the sparsity
  // timing is negligible
  //TIMER tripletsTimer("Building matrix sparsity");
  int DOFs = _vertices.size() * 3;
  _sparseA = SPARSE_MATRIX(DOFs, DOFs);
  _sparseA.setFromTriplets(triplets.begin(), triplets.end());
  _sparseA.makeCompressed();
  //tripletsTimer.stop();
  cout << "done. " << endl;

  // find all the compressed index mapping
  computeCompressedIndices();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::initializeTets()
{
  TIMER functionTimer(__FUNCTION__);
  computeTetIndices();
  computeDmInvs();
  computePFpxs();

  _restTetVolumes = computeTetVolumes(_restVertices);

  _volumeFs.resize(_totalTets);
  _twistFs.resize(_totalTets);
  _abWeights.resize(_totalTets);
  _twistPFPxs.resize(_totalTets);
  _edgeFs.resize(_totalEdges);
  _bendingEs.resize(_totalBends);
  _Us.resize(_totalTets);
  _Vs.resize(_totalTets);
  _Sigmas.resize(_totalTets);

  // is the vertex at the end of a strand?
  _isStrandEnd.clear();
  _isStrandBegin.clear();
  for (unsigned int x = 0; x < _totalVertices; x++)
  {
    _isStrandEnd.push_back(false);
    _isStrandBegin.push_back(false);
  }
  for (unsigned int x = 0; x < _strandIndices.size(); x++)
  {
    const vector<int>& strand = _strandIndices[x];
    _isStrandEnd[strand.back()] = true;
    _isStrandBegin[strand.front()] = true;
  }

  _materials.resize(_totalTets);
  for (int x = 0; x < _totalTets; x++)
  {
    const VECTOR4I tet = _tets[x];

#if 0
    //const REAL mu = 10000.0;  // creates an inversion artifact, if there's only one volume energy per tet
    const REAL mu = 5000.0; // with two volume terms, artifact appears
    //const REAL mu = 1000.0;
    //const REAL mu = 500.0;
    //const REAL mu = 100.0;
    //const REAL mu = 10.0;
    const REAL ks = mu;
    const REAL kb = mu;
    const REAL kt = mu;
#else
    const REAL ks = _E * M_PI * _radiusA * _radiusB;
    const REAL kb = _E * M_PI * (_radiusA * _radiusA * _radiusA * _radiusB) * 0.25;
    //const REAL kt = ((1.0 / 8.0) * _E / (1.0 + _nu)) * _radiusA * _radiusB *
    //                (_radiusA * _radiusA + _radiusB * _radiusB);
    const REAL kt = (M_PI / 4.0) * _G * _radiusA * _radiusB *
                    (_radiusA * _radiusA + _radiusB * _radiusB);
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

    // add twist energy for base triangle (twist energy replacements actually done in TET_WISP_MESH)
    const VECTOR3 normal0 = (v1 - v0).cross(v2 - v0).normalized();
    const VECTOR3 normal1 = (v1 - v2).cross(v3 - v2).normalized();
    _materials[x].volumeEnergies().push_back(new VOLUME::ANISOTROPIC_ARAP(kt, normal0));
    _materials[x].volumeEnergies().push_back(new VOLUME::ANISOTROPIC_ARAP(kt, normal1));
    _materials[x].volumeEnergies().push_back(new VOLUME::ANISOTROPIC_DIRICHLET(kt, normal0));
    _materials[x].volumeEnergies().push_back(new VOLUME::ANISOTROPIC_DIRICHLET(kt, normal1));

    // if it's the first tet in a strand, stretching needs to be handled here too
    if (_isStrandBegin[tet[0]])
    {
      _materials[x].stretchingEnergies()[0] = new STRAND::QUADRATIC_STRETCHING(ks);
      _materials[x].stretchingEnergies()[1] = new STRAND::QUADRATIC_STRETCHING(ks);
    }
    // if it's the last tet in a strand, stretching needs to be handled here too
    if (_isStrandEnd[tet[3]])
    {
      const VECTOR3& v3 = _restVertices[tet[3]];
      MATRIX3x2 tipE;
      tipE.col(0) = v1 - v2;
      tipE.col(1) = v3 - v2;
      const REAL tipTheta0 = STRAND::ISOTROPIC_BENDING::angle(tipE);
      //_materials[x].bendingEnergies()[1] = new STRAND::QUADRATIC_BENDING(kb, tipTheta0);
      _materials[x].bendingEnergies()[1] = new STRAND::QUADRATIC_UNIT_BENDING(kb, tipTheta0);
    }
  }
  computeFs();

  // DEBUG: peek at the forces
  _perTetForces.resize(_totalTets);
}
   
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::computeTetIndices()
{
  // clear any old tets
  _tets.clear();
  _tetEdges.clear();

  // build a map from edge pairs to the corresponding bends
  map<pair<int,int>, int> bendHash;
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    pair<int,int> hash(_bendEdges[x][0], _bendEdges[x][1]);
    bendHash[hash] = x;
  }

  // for each strand
  for (unsigned int y = 0; y < _totalStrands; y++)
  {
    const vector<int>& strand = _strandIndices[y];
    // there's more than 4 points in a strand, right?
    assert(strand.size() >= 4);

    // make each successive 4-tuple a tet
    for (unsigned int x = 3; x < strand.size(); x++)
    {
      VECTOR4I tet;
      tet[0] = strand[x - 3];
      tet[1] = strand[x - 2];
      tet[2] = strand[x - 1];
      tet[3] = strand[x];
      _tets.push_back(tet);

      // which edges are contained inside this tet?
      VECTOR3I edges;
      const int e0 = _edgeHash.find(pair<int,int>(tet[0], tet[1]))->second;
      const int e1 = _edgeHash.find(pair<int,int>(tet[1], tet[2]))->second;
      const int e2 = _edgeHash.find(pair<int,int>(tet[2], tet[3]))->second;
      edges[0] = e0;
      edges[1] = e1;
      edges[2] = e2;
      _tetEdges.push_back(edges);

      VECTOR2I bends;
      bends[0] = bendHash.find(pair<int, int>(e0, e1))->second;
      bends[1] = bendHash.find(pair<int, int>(e1, e2))->second;
      _tetBends.push_back(bends);
    }
  }
  _totalTets = _tets.size();
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::computeDmInvs()
{
  _volumeDmInvs.resize(_totalTets);
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    MATRIX3 Dm;
    Dm.col(0) = _vertices[tet[1]] - _vertices[tet[0]];
    Dm.col(1) = _vertices[tet[2]] - _vertices[tet[0]];
    Dm.col(2) = _vertices[tet[3]] - _vertices[tet[0]];
    _volumeDmInvs[i] = Dm.inverse();
  }

  _edgeDmInvs.resize(_totalEdges);
  for (unsigned int x = 0; x < _totalEdges; x++)
    _edgeDmInvs[x] = 1.0 / _restEdgeLengths[x];
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::computePFpxs()
{
  _volumePFPxs.resize(_totalTets);
  for (unsigned int i = 0; i < _volumeDmInvs.size(); i++)
  {
    const MATRIX3& DmInv = _volumeDmInvs[i];

    const REAL m = DmInv(0, 0);
    const REAL n = DmInv(0, 1);
    const REAL o = DmInv(0, 2);
    const REAL p = DmInv(1, 0);
    const REAL q = DmInv(1, 1);
    const REAL r = DmInv(1, 2);
    const REAL s = DmInv(2, 0);
    const REAL t = DmInv(2, 1);
    const REAL u = DmInv(2, 2);

    const REAL t1 = -m - p - s;
    const REAL t2 = -n - q - t;
    const REAL t3 = -o - r - u;

    MATRIX9x12 PFPu = MATRIX9x12::Zero();
    PFPu(0, 0)  = t1;
    PFPu(0, 3)  = m;
    PFPu(0, 6)  = p;
    PFPu(0, 9)  = s;
    PFPu(1, 1)  = t1;
    PFPu(1, 4)  = m;
    PFPu(1, 7)  = p;
    PFPu(1, 10) = s;
    PFPu(2, 2)  = t1;
    PFPu(2, 5)  = m;
    PFPu(2, 8)  = p;
    PFPu(2, 11) = s;
    PFPu(3, 0)  = t2;
    PFPu(3, 3)  = n;
    PFPu(3, 6)  = q;
    PFPu(3, 9)  = t;
    PFPu(4, 1)  = t2;
    PFPu(4, 4)  = n;
    PFPu(4, 7)  = q;
    PFPu(4, 10) = t;
    PFPu(5, 2)  = t2;
    PFPu(5, 5)  = n;
    PFPu(5, 8)  = q;
    PFPu(5, 11) = t;
    PFPu(6, 0)  = t3;
    PFPu(6, 3)  = o;
    PFPu(6, 6)  = r;
    PFPu(6, 9)  = u;
    PFPu(7, 1)  = t3;
    PFPu(7, 4)  = o;
    PFPu(7, 7)  = r;
    PFPu(7, 10) = u;
    PFPu(8, 2)  = t3;
    PFPu(8, 5)  = o;
    PFPu(8, 8)  = r;
    PFPu(8, 11) = u;

    _volumePFPxs[i] = PFPu;
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::computeTwistPFpxs()
{
  for (unsigned int i = 0; i < _totalTets; i++)
  {
    const REAL alph = _abWeights[i].first;
    const REAL bet = _abWeights[i].second;

    MATRIX9x12 pFpX;
    pFpX.setZero();
    pFpX(0,0) = 1;
    pFpX(0,3) = -1 + alph;
    pFpX(0,6) = -alph;
    pFpX(1,1) = 1;
    pFpX(1,4) = -1 + alph;
    pFpX(1,7) = -alph;
    pFpX(2,2) = 1;
    pFpX(2,5) = -1 + alph;
    pFpX(2,8) = -alph;
    pFpX(3,3) = bet;
    pFpX(3,6) = -1 - bet;
    pFpX(3,9) = 1;
    pFpX(4,4) = bet;
    pFpX(4,7) = -1 - bet;
    pFpX(4,10) = 1;
    pFpX(5,5) = bet;
    pFpX(5,8) = -1 - bet;
    pFpX(5,11) = 1;
    pFpX(6,3) = -1;
    pFpX(6,6) = 1;
    pFpX(7,4) = -1;
    pFpX(7,7) = 1;
    pFpX(8,5) = -1;
    pFpX(8,8) = 1;

    _twistPFPxs[i] = pFpX;
  }
}

///////////////////////////////////////////////////////////////////////
// compute volumes for tets -- works for rest and deformed, just
// pass it _restVertices or _vertices
///////////////////////////////////////////////////////////////////////
VECTOR TET_STRAND_MESH::computeTetVolumes(const vector<VECTOR3>& vertices)
{
  VECTOR tetVolumes(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    vector<VECTOR3> tetVertices;
    for (unsigned int j = 0; j < 4; j++)
      tetVertices.push_back(_vertices[tet[j]]);
    tetVolumes[i] = computeTetVolume(tetVertices);

    if (tetVolumes[i] < 0.0)
    {
      cout << " ERROR: Bad rest volume found: " << tetVolumes[i] << endl;
      cout << " Tet index: " << i << endl;
      // tetVolumes[i] = 0.0;
    }
    assert(tetVolumes[i] >= 0.0);
  }
  return tetVolumes;
}

///////////////////////////////////////////////////////////////////////
// compute the volume of a tet
///////////////////////////////////////////////////////////////////////
REAL TET_STRAND_MESH::computeTetVolume(const vector<VECTOR3>& tetVertices)
{
  const VECTOR3 diff1 = tetVertices[1] - tetVertices[0];
  const VECTOR3 diff2 = tetVertices[2] - tetVertices[0];
  const VECTOR3 diff3 = tetVertices[3] - tetVertices[0];
  return diff3.dot((diff1).cross(diff2)) / 6.0;
}

// change-of-basis for edges on a tet to spatial gradient
static std::array<MATRIX3x12, 3> edgeDfDx()
{
  std::array<MATRIX3x12,3> result;
  result[0].setZero();
  result[1].setZero();
  result[2].setZero();

  MATRIX3 I = MATRIX3::Identity();
  result[0].block<3,3>(0,0) = -I;
  result[0].block<3,3>(0,3) = I;
  
  result[1].block<3,3>(0,3) = -I;
  result[1].block<3,3>(0,6) = I;
  
  result[2].block<3,3>(0,6) = -I;
  result[2].block<3,3>(0,9) = I;

  return result;
}

// change-of-basis for edges on a tet to spatial gradient
static std::array<MATRIX, 2> bendDfDx()
{
  std::array<MATRIX,2> result;
  result[0] = MATRIX(6,12);
  result[1] = MATRIX(6,12);
  result[0].setZero();
  result[1].setZero();

  MATRIX3 I = MATRIX3::Identity();
  result[0].block<3,3>(0,0) =  I;
  result[0].block<3,3>(0,3) = -I;
  result[0].block<3,3>(3,3) = -I;
  result[0].block<3,3>(3,6) =  I;
  
  result[1].block<3,3>(0,3) =  I;
  result[1].block<3,3>(0,6) = -I;
  result[1].block<3,3>(3,6) = -I;
  result[1].block<3,3>(3,9) =  I;
  return result;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
REAL TET_STRAND_MESH::computeStretchingEnergy()
{
  REAL energy = 0;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR3I edges = _tetEdges[i];
    std::array<VECTOR3,3> fs;

    fs[0] = _edgeFs[edges[0]];
    fs[1] = _edgeFs[edges[1]];
    fs[2] = _edgeFs[edges[2]];
    
    VECTOR3 restLengths;
    restLengths[0] = _restEdgeLengths[edges[0]];
    restLengths[1] = _restEdgeLengths[edges[1]];
    restLengths[2] = _restEdgeLengths[edges[2]];

    const std::array<REAL,3> energies = _materials[i].edgeEnergies(fs);
    for (unsigned int j = 0; j < 3; j++)
      energy += restLengths[j] * energies[j];
  }
  
  return energy;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TET_STRAND_MESH::computeStretchingForces()
{
  //TIMER functionTimer(__FUNCTION__);

  // get the edge forces
  vector<VECTOR12> perElementForces(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
    perElementForces[i].setZero();
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR3I edges = _tetEdges[i];
    std::array<VECTOR3,3> fs;

    fs[0] = _edgeFs[edges[0]];
    fs[1] = _edgeFs[edges[1]];
    fs[2] = _edgeFs[edges[2]];
    
    VECTOR3 restLengths;
    restLengths[0] = _restEdgeLengths[edges[0]];
    restLengths[1] = _restEdgeLengths[edges[1]];
    restLengths[2] = _restEdgeLengths[edges[2]];

    const std::array<VECTOR3,3> PK1s = _materials[i].edgePK1(fs);
    const std::array<MATRIX3x12, 3> dfdxs = edgeDfDx();

    for (unsigned int j = 0; j < 3; j++)
    {
      const VECTOR12 force = -restLengths[j] * (dfdxs[j].transpose() * PK1s[j]);
      perElementForces[i] += force;
    }
    
    if (perElementForces[i].hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " Found NaN in tet " << i << endl;
      cout << " PK1: " << endl << PK1s[0] << endl;
      cout << " PK1: " << endl << PK1s[1] << endl;
      cout << " PK1: " << endl << PK1s[2] << endl;
      cout << " edges: " << endl << fs[0].transpose() 
                         << endl << fs[1].transpose() 
                         << endl << fs[2].transpose() << endl;
      cout << " vertices: " << _tets[i].transpose() << endl;
      exit(0);
    }
  }
  
  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int y = 0; y < _tets.size(); y++)
  {
    const VECTOR4I& tet = _tets[y];
    const VECTOR12& tetForce = perElementForces[y];
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
///////////////////////////////////////////////////////////////////////
REAL TET_STRAND_MESH::computeBendingEnergy()
{
  REAL energy = 0;

  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR2I bend = _tetBends[i];
    std::array<MATRIX3x2,2> Es;
    Es[0] = _bendingEs[bend[0]];
    Es[1] = _bendingEs[bend[1]];

    const std::array<REAL,2> energies = _materials[i].bendEnergies(Es);
    for (int j = 0; j < 2; j++)
    {
      const REAL area = _voronoiLengths[bend[j]];
      energy += area * energies[j];
    }
  }

  return energy;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TET_STRAND_MESH::computeBendingForces()
{
  //TIMER functionTimer(__FUNCTION__);

  // get the edge forces
  vector<VECTOR12> perElementForces(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
    perElementForces[i].setZero();
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR2I bend = _tetBends[i];
    std::array<MATRIX3x2,2> Es;
    Es[0] = _bendingEs[bend[0]];
    Es[1] = _bendingEs[bend[1]];

    const std::array<MATRIX3x2,2> PK1s = _materials[i].bendPK1(Es);
    const std::array<MATRIX, 2> dfdxs = bendDfDx();

    for (int j = 0; j < 2; j++)
    {
      const REAL area = _voronoiLengths[bend[j]];
      const VECTOR12 force = -area * dfdxs[j].transpose() * flatten(PK1s[j]);
      perElementForces[i] += force;
    }
    
    if (perElementForces[i].hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " Found NaN in tet " << i << endl;
      cout << " Bend 0: " << endl << Es[0] << endl;
      cout << " Bend 1: " << endl << Es[1] << endl;
      cout << " PK1 0: " << endl << PK1s[0] << endl;
      cout << " PK1 1: " << endl << PK1s[1] << endl;
      exit(0);
    }
  }

  // DEBUG
  _perTetForces = perElementForces;

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int y = 0; y < _tets.size(); y++)
  {
    const VECTOR4I& tet = _tets[y];
    const VECTOR12& tetForce = perElementForces[y];
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
///////////////////////////////////////////////////////////////////////
REAL TET_STRAND_MESH::computeTwistingEnergy()
{
  REAL energy = 0;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const MATRIX3& twistF = _twistFs[i];
    const REAL twistDensity = _materials[i].twistEnergy(twistF);
    const REAL twistEnergy = _restTetVolumes[i] * twistDensity;
    //cout << " twist density " << i << ": " << twistDensity<< endl;
    //cout << " rest volume: " << _restTetVolumes[i] << endl;
    //cout << " twist energy " << i << ": " << twistEnergy << endl;
    energy += twistEnergy;
  }
  
  return energy;
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
VECTOR TET_STRAND_MESH::computeTwistingForces()
{
  //TIMER functionTimer(__FUNCTION__);

  // get the edge forces
  vector<VECTOR12> perElementForces(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
    perElementForces[i].setZero();
  
  computeTwistPFpxs();

//#pragma omp parallel
//#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const MATRIX3& twistF = _twistFs[i];
    const MATRIX3& volumeF = _volumeFs[i];
#if 0
    if (F.determinant() < 0.0)
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " INVERTED TET: " << i << endl;
      cout << " F: " << endl << F << endl;
    }
#endif
    const MATRIX3 volumePK1 = _materials[i].volumePK1(volumeF);
    const MATRIX3 twistPK1 = _materials[i].twistPK1(twistF);

    const VECTOR12 volumeDensity = _volumePFPxs[i].transpose() * flatten(volumePK1);
    const VECTOR12 twistDensity = _twistPFPxs[i].transpose() * flatten(twistPK1);
    const VECTOR12 force = -_restTetVolumes[i] * (volumeDensity + twistDensity);
    if (force.hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " Found NaN in tet " << i << endl;
      cout << " volumeF: " << endl << volumeF << endl;
      cout << " volumePK1: " << endl << volumePK1 << endl;
      cout << " twistF: " << endl << twistF << endl;
      cout << " twistPK1: " << endl << twistPK1 << endl;
    }
    perElementForces[i] = force;
  }
  // DEBUG
  //_perTetForces = perElementForces;

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int y = 0; y < _tets.size(); y++)
  {
    const VECTOR4I& tet = _tets[y];
    const VECTOR12& tetForce = perElementForces[y];
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
///////////////////////////////////////////////////////////////////////
VECTOR TET_STRAND_MESH::computeHyperelasticForces()
{
  TIMER functionTimer(__FUNCTION__);

  // get the volume forces
  vector<VECTOR12> perElementForces(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    // TOGGLE HERE for edge twisting 
    // const MATRIX3& F = _volumeFs[i];
    const MATRIX3& F = _twistFs[i];
    if (F.determinant() < 0.0)
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " INVERTED TET: " << i << endl;
    }
    const MATRIX3 PK1 = _materials[i].volumePK1(F);
    // TOGGLE here for edge twisting
    // const VECTOR12 forceDensity = _volumePFPxs[i].transpose() * flatten(PK1);
    const VECTOR12 forceDensity = _twistPFPxs[i].transpose() * flatten(PK1);
    const VECTOR12 force = -_restTetVolumes[i] * forceDensity;
    if (force.hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " Found NaN in tet " << i << endl;
      cout << " F: " << endl << F << endl;
      cout << " PK1: " << endl << PK1 << endl;
    }
    perElementForces[i] = force;
  }

  // DEBUG, to peek at forces
  _perTetForces = perElementForces;

  // get the edge forces
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR3I edges = _tetEdges[i];
    std::array<VECTOR3,3> fs;

    fs[0] = _edgeFs[edges[0]];
    fs[1] = _edgeFs[edges[1]];
    fs[2] = _edgeFs[edges[2]];
    
    VECTOR3 restLengths;
    restLengths[0] = _restEdgeLengths[edges[0]];
    restLengths[1] = _restEdgeLengths[edges[1]];
    restLengths[2] = _restEdgeLengths[edges[2]];

    const std::array<VECTOR3,3> PK1s = _materials[i].edgePK1(fs);
    const std::array<MATRIX3x12, 3> dfdxs = edgeDfDx();

    for (unsigned int j = 0; j < 3; j++)
    {
      const VECTOR12 force = -restLengths[j] * (dfdxs[j].transpose() * PK1s[j]);
      perElementForces[i] += force;
    }
    
    if (perElementForces[i].hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " Found NaN in tet " << i << endl;
    }
  }

  // get the bend forces
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR2I bend = _tetBends[i];
    std::array<MATRIX3x2,2> Es;
    Es[0] = _bendingEs[bend[0]];
    Es[1] = _bendingEs[bend[1]];

    const std::array<MATRIX3x2,2> PK1s = _materials[i].bendPK1(Es);
    const std::array<MATRIX, 2> dfdxs = bendDfDx();

    for (int j = 0; j < 2; j++)
    {
      const REAL area = _voronoiLengths[bend[j]];
      const VECTOR12 force = -area * dfdxs[j].transpose() * flatten(PK1s[j]);
      perElementForces[i] += force;
    }
    
    if (perElementForces[i].hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " Found NaN in tet " << i << endl;
    }
  }

  // DEBUG, to peek at forces
  //_perTetForces = perElementForces;

  // scatter the forces to the global force vector, this can be parallelized
  // better where each vector entry pulls from perElementForce, but let's get
  // the slow preliminary version working first
  const int DOFs = _vertices.size() * 3;
  VECTOR forces(DOFs);
  forces.setZero();

  for (unsigned int y = 0; y < _tets.size(); y++)
  {
    const VECTOR4I& tet = _tets[y];
    const VECTOR12& tetForce = perElementForces[y];
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
// redo fiber direction for each tet
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::recomputeTwistNormals()
{
  for (int x = 0; x < _totalTets; x++)
  {
    const VECTOR4I tet = _tets[x];
    // add the bending energy 
    const VECTOR3& v0 = _restVertices[tet[0]];
    const VECTOR3& v1 = _restVertices[tet[1]];
    const VECTOR3& v2 = _restVertices[tet[2]];
    const VECTOR3& v3 = _restVertices[tet[3]];

    //const REAL kt = ((1.0 / 8.0) * _E / (1.0 + _nu)) * _radiusA * _radiusB *
    //                (_radiusA * _radiusA + _radiusB * _radiusB);
    const REAL kt = (M_PI / 4.0) * _G * _radiusA * _radiusB *
                    (_radiusA * _radiusA + _radiusB * _radiusB);

    // add twist energy for base triangle
    const VECTOR3 normal0 = (v1 - v0).cross(v2 - v0).normalized();
    const VECTOR3 normal1 = (v1 - v2).cross(v3 - v2).normalized();
    const MATRIX3& F = _volumeFs[x];
    const VECTOR3 nHat0 = (F.inverse() * normal0).normalized();
    const VECTOR3 nHat1 = (F.inverse() * normal1).normalized();
    _materials[x].volumeEnergies().clear();
    _materials[x].volumeEnergies().push_back(new VOLUME::ANISOTROPIC_ARAP(kt, nHat0));
    _materials[x].volumeEnergies().push_back(new VOLUME::ANISOTROPIC_ARAP(kt, nHat1));
    _materials[x].volumeEnergies().push_back(new VOLUME::ANISOTROPIC_DIRICHLET(kt, nHat0));
    _materials[x].volumeEnergies().push_back(new VOLUME::ANISOTROPIC_DIRICHLET(kt, nHat1));
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::computeFs()
{
  //TIMER functionTimer(__FUNCTION__);
  assert(_volumeFs.size() == (unsigned int)_totalTets);
  for (int x = 0; x < _totalTets; x++)
  {
    const VECTOR4I& tet = _tets[x];
    MATRIX3 Ds;
    Ds.col(0) = _vertices[tet[1]] - _vertices[tet[0]];
    Ds.col(1) = _vertices[tet[2]] - _vertices[tet[0]];
    Ds.col(2) = _vertices[tet[3]] - _vertices[tet[0]];
    _volumeFs[x] = Ds * _volumeDmInvs[x];

    // Filling in _twistFs and the _abWeights
    VECTOR3 e1 = _vertices[tet[2]] - _vertices[tet[1]];
    VECTOR3 e0 = _vertices[tet[0]] - _vertices[tet[1]];
    VECTOR3 e2 = _vertices[tet[3]] - _vertices[tet[2]];

    _abWeights[x].first = e0.dot(e1)/e1.squaredNorm();
    _abWeights[x].second = e2.dot(e1)/e1.squaredNorm();
    
    e0 = e0 - _abWeights[x].first * e1;
    e2 = e2 - _abWeights[x].second * e1;
    MATRIX3 Ftwist;
    Ftwist.col(0) = e0;
    Ftwist.col(1) = e2;
    Ftwist.col(2) = e1;
    _twistFs[x] = Ftwist;
  }
  
  for (unsigned int x = 0; x < _totalEdges; x++)
  {
    const VECTOR2I edge = _edgeIndices[x];
    const VECTOR3 v0 = _vertices[edge[0]];
    const VECTOR3 v1 = _vertices[edge[1]];

    _edgeFs[x] = (v1 - v0) * _edgeDmInvs[x];
  }
  
  for (unsigned int x = 0; x < _totalBends; x++)
  {
    const VECTOR3I bend = _bendVertices[x];
    const VECTOR3& v0 = _vertices[bend[0]];
    const VECTOR3& v1 = _vertices[bend[1]];
    const VECTOR3& v2 = _vertices[bend[2]];

    MATRIX3x2 E;
    E.col(0) = v0 - v1;
    E.col(1) = v2 - v1;
    _bendingEs[x] = E;
  }
  
  // record which ones are inverted
  _isTetInverted.resize(_totalTets);
  for (int x = 0; x < _totalTets; x++)
    _isTetInverted[x] = _volumeFs[x].determinant() < 0.0;

  // recompute edge length for inversion handling
  computeEdgeLengths();

  // DEBUG: recompute normal directions
  //recomputeTwistNormals();
}

///////////////////////////////////////////////////////////////////////
// get the SVD of all the deformation gradients
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::computeSVDs()
{
  //TIMER functionTimer(__FUNCTION__);
  assert(_Us.size() == _tets.size());
  assert(_Sigmas.size() == _tets.size());
  assert(_Vs.size() == _tets.size());
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < _tets.size(); x++)
    svd_rv(_volumeFs[x], _Us[x], _Sigmas[x], _Vs[x]);
}

///////////////////////////////////////////////////////////////////////
// scale strand material parameters
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::scaleStretchingMu(const REAL& scale)
{
  for (unsigned int x = 0; x < _materials.size(); x++)
  {
    using namespace std;
    using namespace STRAND;
    array<STRETCHING*,3>& stretchingEnergies = _materials[x].stretchingEnergies();

    for (unsigned int y = 0; y < 3; y++)
    {
      if (stretchingEnergies[y] == NULL) continue;
      stretchingEnergies[y]->mu() *= scale;
    }
  }
}

///////////////////////////////////////////////////////////////////////
// scale strand material parameters
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::scaleBendingMu(const REAL& scale)
{
  for (unsigned int x = 0; x < _materials.size(); x++)
  {
    using namespace std;
    using namespace STRAND;
    array<ISOTROPIC_BENDING*,2>& bendingEnergies = _materials[x].bendingEnergies();

    for (unsigned int y = 0; y < 2; y++)
    {
      if (bendingEnergies[y] == NULL) continue;
      bendingEnergies[y]->mu() *= scale;
    }
  }
}

///////////////////////////////////////////////////////////////////////
// scale strand material parameters
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::scaleTwistingMu(const REAL& scale)
{
  for (unsigned int x = 0; x < _materials.size(); x++)
  {
    vector<VOLUME::HYPERELASTIC*>& volumeEnergies = _materials[x].volumeEnergies();

    for (unsigned int y = 0; y < volumeEnergies.size(); y++)
    {
      using namespace VOLUME;
      TET_STRAND_TWIST* twist = (TET_STRAND_TWIST*)volumeEnergies[y];
      twist->mu() *= scale;
    }
  }
}

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_STRAND_MESH::computeHyperelasticClampedHessian()
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perElementHessians(_tets.size());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    //const MATRIX3& twistF       = _twistFs[i];
    //const MATRIX9x12& twistPFpx = _twistPFPxs[i];
    const MATRIX3& volumeF       = _volumeFs[i];
    const MATRIX9x12& volumePFpx = _volumePFPxs[i];
    const MATRIX9 volumeHessian = -_restTetVolumes[i] * _materials[i].volumeClampedHessian(volumeF);
    //const MATRIX9 twistHessian = -_restTetVolumes[i] * _materials[i].twistClampedHessian(twistF);

    perElementHessians[i].setZero();
    vector<VECTOR3> vertices(4);
    vertices[0] = _vertices[_tets[i][0]];
    vertices[1] = _vertices[_tets[i][1]];
    vertices[2] = _vertices[_tets[i][2]];
    vertices[3] = _vertices[_tets[i][3]];

    const MATRIX12 twistHessian = -_restTetVolumes[i] * _materials[i].twistClampedForceGradient(vertices);
    perElementHessians[i] += twistHessian;
    perElementHessians[i] += (volumePFpx.transpose() * volumeHessian) * volumePFpx;
    //perElementHessians[i] += (twistPFpx.transpose() * twistHessian) * twistPFpx;
   
    // edge terms 
    const VECTOR3I edges = _tetEdges[i];
    std::array<VECTOR3,3> fs;

    fs[0] = _edgeFs[edges[0]];
    fs[1] = _edgeFs[edges[1]];
    fs[2] = _edgeFs[edges[2]];
    
    VECTOR3 restLengths;
    restLengths[0] = _restEdgeLengths[edges[0]];
    restLengths[1] = _restEdgeLengths[edges[1]];
    restLengths[2] = _restEdgeLengths[edges[2]];

    const std::array<MATRIX3,3> edgeHs = _materials[i].edgeClampedHessian(fs);
    const std::array<MATRIX3x12, 3> edgeDfdxs = edgeDfDx();

    for (unsigned int j = 0; j < 3; j++)
      perElementHessians[i] += -restLengths[j] * (edgeDfdxs[j].transpose() * edgeHs[j] * edgeDfdxs[j]);
   
    // bend terms 
    const VECTOR2I bend = _tetBends[i];
    std::array<MATRIX3x2,2> Es;
    Es[0] = _bendingEs[bend[0]];
    Es[1] = _bendingEs[bend[1]];
    const std::array<MATRIX6,2> bendHs = _materials[i].bendClampedHessian(Es);
    const std::array<MATRIX, 2> bendDfdxs = bendDfDx();

    for (int j = 0; j < 2; j++)
    {
      const REAL area = _voronoiLengths[bend[j]];
      perElementHessians[i] += -area * bendDfdxs[j].transpose() * bendHs[j] * bendDfdxs[j];
    }
  }

  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;
  vector<TRIPLET> triplets;
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
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
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_STRAND_MESH::computeHyperelasticClampedHessianFast()
{
  TIMER functionTimer(__FUNCTION__);
  vector<MATRIX12> perElementHessians(_tets.size());
//cout << " OMP DISABLED " << __FUNCTION__ << endl; 
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const MATRIX3& twistF       = _twistFs[i];
    const MATRIX9x12& twistPFpx = _twistPFPxs[i];
    const MATRIX3& volumeF       = _volumeFs[i];
    const MATRIX9x12& volumePFpx = _volumePFPxs[i];
    const MATRIX9 volumeHessian = -_restTetVolumes[i] * _materials[i].volumeClampedHessian(volumeF);
    perElementHessians[i] = (volumePFpx.transpose() * volumeHessian) * volumePFpx;

#if 0
    const MATRIX9 twistHessian = -_restTetVolumes[i] * _materials[i].twistClampedHessian(twistF);
    perElementHessians[i] += (twistPFpx.transpose() * twistHessian) * twistPFpx;
#else
    vector<VECTOR3> vertices(4);
    vertices[0] = _vertices[_tets[i][0]];
    vertices[1] = _vertices[_tets[i][1]];
    vertices[2] = _vertices[_tets[i][2]];
    vertices[3] = _vertices[_tets[i][3]];

    const MATRIX12 twistHessian = -_restTetVolumes[i] * _materials[i].twistClampedForceGradient(vertices);
    perElementHessians[i] += twistHessian;
#endif

    if (perElementHessians[i].hasNaN())
    {
      std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
      cout << " Found NaN in tet " << i << " volume term " << endl;

      if (volumeHessian.hasNaN())
        cout << " volume term has NaN " << endl;
      if (twistHessian.hasNaN())
      {
        cout << " twist term has NaN " << endl;
        cout << " twist F: " << endl << twistF << endl;
        cout << " twist hessian: " << endl << twistHessian << endl;
      }
      
      exit(0);
    }

    // edge terms 
    const VECTOR3I edges = _tetEdges[i];
    std::array<VECTOR3,3> fs;

    fs[0] = _edgeFs[edges[0]];
    fs[1] = _edgeFs[edges[1]];
    fs[2] = _edgeFs[edges[2]];
    
    VECTOR3 restLengths;
    restLengths[0] = _restEdgeLengths[edges[0]];
    restLengths[1] = _restEdgeLengths[edges[1]];
    restLengths[2] = _restEdgeLengths[edges[2]];

    const std::array<MATRIX3,3> edgeHs = _materials[i].edgeClampedHessian(fs);
    const std::array<MATRIX3x12, 3> edgeDfdxs = edgeDfDx();

    for (unsigned int j = 0; j < 3; j++)
    {
      perElementHessians[i] += -restLengths[j] * (edgeDfdxs[j].transpose() * edgeHs[j] * edgeDfdxs[j]);
      if (edgeHs[j].hasNaN())
      {
        std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
        cout << " Found NaN in tet " << i << " edge term " << endl;
        exit(0);
      }
    }
    

    // bend terms 
    const VECTOR2I bend = _tetBends[i];
    std::array<MATRIX3x2,2> Es;
    Es[0] = _bendingEs[bend[0]];
    Es[1] = _bendingEs[bend[1]];
    const std::array<MATRIX6,2> bendHs = _materials[i].bendClampedHessian(Es);
    const std::array<MATRIX, 2> bendDfdxs = bendDfDx();

    for (int j = 0; j < 2; j++)
    {
      const REAL area = _voronoiLengths[bend[j]];
      perElementHessians[i] += -area * bendDfdxs[j].transpose() * bendHs[j] * bendDfdxs[j];
      if (bendHs[j].hasNaN())
      {
        std::cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << std::endl;
        cout << " Found NaN in tet " << i << " bend term " << endl;
        cout << " Es: " << endl << Es[j] << endl;
        cout << " bendHs: " << endl << bendHs[j] << endl;
        exit(0);
      }
    }
  }

  // could probably do better here by:
  // 1. arranging things into 3x3 blocks instead of entry-wise
  // 2. using symmetry so we don't set the same entry twice
  // this isn't at the top of the timing pile anymore though.
  TIMER assemblyTimer("Sparse matrix assembly");
  const unsigned int nonZeros = _sparseA.nonZeros();
  REAL* base = _sparseA.valuePtr();
#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int x = 0; x < nonZeros; x++)
  {
    const vector<VECTOR3I>& gather = _hessianGathers[x];
    base[x] = 0;

    for (unsigned int y = 0; y < gather.size(); y++)
    {
      const VECTOR3I& lookup = gather[y];
      const int& tetIndex = lookup[0];
      const int& row = lookup[1];
      const int& col = lookup[2];
      base[x] += perElementHessians[tetIndex](row, col);
    }
  }
  return _sparseA;
}

///////////////////////////////////////////////////////////////////////
// set the vertex positions directly exactly
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::setDisplacement(const VECTOR& displacements)
{
  assert((unsigned int)displacements.size() == _DOFs);

  // back up old positions, for visualization
  _verticesOld = _vertices;

  unsigned int index = 0;
  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    for (unsigned int y = 0; y < 3; y++, index++)
      _vertices[x][y] = _restVertices[x][y] + displacements[index];
  }
}

///////////////////////////////////////////////////////////////////////
// get the vertex positions, flattened into a vector
///////////////////////////////////////////////////////////////////////
const VECTOR TET_STRAND_MESH::getDisplacement() const
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

///////////////////////////////////////////////////////////////////////
// find the compressed index mapping
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::computeCompressedIndices()
{
  TIMER functionTimer(__FUNCTION__);

  cout << " Hashing indices ... " << flush;

  // cache the beginning of the storage
  REAL* base = _sparseA.valuePtr();

  for (unsigned int x = 0; x < _sparseA.outerSize(); x++)
  {
    for (SPARSE_MATRIX::InnerIterator it(_sparseA, x); it; ++it)
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

  cout << " Computing compressed indices ... " << flush;
  // allocate an array for each non-zero matrix entry
  _hessianGathers.resize(_sparseA.nonZeros());
  for (unsigned int i = 0; i < _tets.size(); i++)
  {
    const VECTOR4I& tet = _tets[i];
    for (int y = 0; y < 4; y++)
    {
      const int yVertex = tet[y];
      for (int x = 0; x < 4; x++)
      {
        const int xVertex = tet[x];
        for (int b = 0; b < 3; b++)
          for (int a = 0; a < 3; a++)
          {
            // do the lookup, see where this is stored globally
            const pair<int, int> rowCol(3 * xVertex + a, 3 * yVertex + b);
            const auto iter = _compressedIndex.find(rowCol);
            const int index = iter->second;

            // store the tet and entry and H this corresponds to
            VECTOR3I tetMapping;
            tetMapping[0] = i;
            tetMapping[1] = 3 * x + a;
            tetMapping[2] = 3 * y + b;
            _hessianGathers[index].push_back(tetMapping);
          }
      }
    }
  }
  cout << "done." << endl;
}

///////////////////////////////////////////////////////////////////////
// edge-edge collision detection
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::computeEdgeEdgeCollisions(const REAL& collisionEps, const bool verbose)
{
#if 0
  STRAND_MESH::computeEdgeEdgeCollisions();
  return;
#endif

  TIMER functionTimer(string("TET_STRAND_MESH::") + string(__FUNCTION__));
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
    const int outerStrand = _perEdgeStrandIndex[x];
    const VECTOR3& v0 = _vertices[outerEdge[0]];
    const VECTOR3& v1 = _vertices[outerEdge[1]];
    //const unsigned int outerFlat = outerEdge[0] + outerEdge[1] * _edgeIndices.size();

    vector<int> nearbyEdges;
    //_collisionTree->nearbyEdges(_edgeIndices[x], _collisionEps, nearbyEdges);
    _collisionTree->nearbyEdges(_edgeIndices[x], collisionEps, nearbyEdges);

    // find the closest other edge
    for (unsigned int y = 0; y < nearbyEdges.size(); y++)
    {
      // skip self
      if ((int)x == nearbyEdges[y]) continue;

      // skip if index is smaller -- don't want to double count nearby edges
      // (a,b) and (b,a)
      if ((unsigned int)nearbyEdges[y] < x) continue;

      const VECTOR2I innerEdge = _edgeIndices[nearbyEdges[y]];
      const int innerStrand = _perEdgeStrandIndex[nearbyEdges[y]];

      // if we're ignoring intra-strand collisions, check here
      if (_strandSelfCollisionDisabled && (innerStrand == outerStrand))
        continue;

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
  if (verbose)
    cout << " Found " << _edgeEdgeCollisions.size() << " edge-edge collisions " << endl;

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
// strand ends are different, since we dropped the twist variable
///////////////////////////////////////////////////////////////////////
void TET_STRAND_MESH::rebuildGlobalIndices()
{
  int index = 0;
  for (unsigned int x = 0; x < _strandIndices.size(); x++)
  {
    const int strandSize = _strandIndices[x].size();
    index += 3 * strandSize;
    _globalStrandEnds[x] = index;
  }
}

///////////////////////////////////////////////////////////////////////
// let's build it all in parallel
///////////////////////////////////////////////////////////////////////
SPARSE_MATRIX TET_STRAND_MESH::buildEdgeEdgeMatrix(const vector<MATRIX12>& perEdgeHessians) const
#if 1
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
// attempted multi-threaded assembly here, but the low-hanging fruit versions don't work.
// Need to go in and build each sparse column in a different thread, which will be more involved.i
// For example, design pattern from: https://stackoverflow.com/questions/47093836/how-to-efficiently-assemble-a-fem-sparse-matrix
#else
{
  TIMER functionTimer(__FUNCTION__);
  // build out the triplets
  typedef Eigen::Triplet<REAL> TRIPLET;

  const int threads = omp_get_max_threads();

  vector<vector<TRIPLET> > triplets(threads);
  for (int x = 0; x < threads; x++)
    triplets[x].reserve(_edgeEdgeCollisions.size());

#pragma omp parallel
#pragma omp for schedule(static)
  for (unsigned int i = 0; i < _edgeEdgeCollisions.size(); i++)
  {
    const MATRIX12& H = perEdgeHessians[i];
    const VECTOR2I& edge0 = _edgeIndices[_edgeEdgeCollisions[i].first];
    const VECTOR2I& edge1 = _edgeIndices[_edgeEdgeCollisions[i].second];

    const int threadID = omp_get_thread_num();

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
            triplets[threadID].push_back(triplet);
          }
      }
    }
  }

#if 1
  SPARSE_MATRIX A(_DOFs, _DOFs);
  vector<SPARSE_MATRIX> As(threads);
#pragma omp parallel
#pragma omp for schedule(static)
  for (int x = 0; x < threads; x++)
  {
    As[x] = SPARSE_MATRIX(_DOFs, _DOFs);
    As[x].setFromTriplets(triplets[x].begin(), triplets[x].end());
#pragma omp critical
    A += As[x];
  }
#else

#pragma omp declare reduction (+: SPARSE_MATRIX: omp_out=omp_out+omp_in)\
    initializer(omp_priv=SPARSE_MATRIX(omp_orig.rows(), omp_orig.cols()))

  TIMER finalAdd("buildEdgeEdgeMatrix: reduction");
#pragma omp parallel for reduction(+: A)
  for (int x = 0; x < threads; x++)
    A += As[x];
  finalAdd.stop();
#endif
  
  //A.setFromTriplets(triplets.begin(), triplets.end());

  return A;
}
#endif

}
