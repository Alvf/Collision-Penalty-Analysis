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
#ifndef STRAND_MESH_H
#define STRAND_MESH_H

#include "SETTINGS.h"
#include "Hyperelastic/Strand/STRETCHING.h"
#include "Hyperelastic/Strand/ISOTROPIC_BENDING.h"
#include "Hyperelastic/Strand/BENDING.h"
//#include "Hyperelastic/Volume/EDGE_COLLISION.h"
#include "Collision/ENERGY_12D.h"

#include <map>
#include <vector>

namespace HOBAK {

using namespace std;

class STRAND_MESH
{
public:
  STRAND_MESH() = default;

  // assumes it's all just one big strand
  STRAND_MESH(const vector<VECTOR3>& restVertices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB);

  // assumes it's all just one big strand
  STRAND_MESH(const vector<VECTOR3>& restVertices,
              const vector<VECTOR3>& vertices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB);

  // accepts a vector of individual strands
  STRAND_MESH(const vector<VECTOR3>& restVertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB);

  // accepts a vector of individual strands, and a modified initial pose
  STRAND_MESH(const vector<VECTOR3>& restVertices,
              const vector<VECTOR3>& vertices,
              const vector<vector<int> >& strandIndices,
              const REAL& E,        // Young's modulus
              const REAL& G,       // shear modulus
              const REAL& density,
              const REAL& radiusA,
              const REAL& radiusB);

  virtual ~STRAND_MESH();

  // accessors
  const vector<VECTOR3>& vertices() const     { return _vertices; };
  const vector<VECTOR3>& verticesOld() const  { return _verticesOld; };
  vector<VECTOR3>& vertices()                 { return _vertices; };
  const vector<VECTOR3>& restVertices() const { return _restVertices; };
  vector<VECTOR3>& restVertices()             { return _restVertices; };
  const vector<VECTOR3>& director1s() const    { return _director1s; };
  vector<VECTOR3>& director1s()                { return _director1s; };
  const vector<VECTOR3>& director2s() const    { return _director2s; };
  vector<VECTOR3>& director2s()                { return _director2s; };
  const vector<VECTOR2I>& edgeIndices() const       { return _edgeIndices; };
  vector<VECTOR2I>& edgeIndices()                   { return _edgeIndices; };
  
  const VECTOR3& vertex(const int index) const     { return _vertices[index]; };
  VECTOR3& vertex(const int index)                 { return _vertices[index]; };
  const VECTOR3& restVertex(const int index) const { return _restVertices[index]; };
  VECTOR3& restVertex(const int index)             { return _restVertices[index]; };
  const VECTOR& restEdgeLengths() const            { return _restEdgeLengths; };
  VECTOR& restEdgeLengths()                        { return _restEdgeLengths; };
  const VECTOR& restTwistAngles() const            { return _restTwistAngles; };
  VECTOR& restTwistAngles()                        { return _restTwistAngles; };
  const VECTOR& restBendAngles() const             { return _restBendAngles; };
  VECTOR& restBendAngles()                         { return _restBendAngles; };
  const VECTOR& referenceTwists() const            { return _referenceTwists; };
  VECTOR& referenceTwists()                        { return _referenceTwists; };
  const VECTOR& twists() const                     { return _twists; };
  VECTOR& twists()                                 { return _twists; };

  const vector<vector<int> >& strandIndices() const { return _strandIndices; };
  vector<vector<int> >& strandIndices()             { return _strandIndices; };
  const unsigned int totalStrands() const           { return _strandIndices.size(); };
  const bool edgeEnd() const                        { return _edgeEnd; };
  const VECTORI& globalVertexIndices() const        { return _globalVertexIndices; };
  const VECTORI& globalEdgeIndices() const          { return _globalEdgeIndices; };
  const VECTORI& globalStrandEnds() const           { return _globalStrandEnds; };

  const vector<pair<int,int> >& edgeEdgeCollisions() const             { return _edgeEdgeCollisions; };
  const vector<pair<VECTOR2,VECTOR2> >& edgeEdgeCoordinates() const    { return _edgeEdgeCoordinates; };
  const vector<pair<int,int> >& edgeEdgeCollisionsOld() const          { return _edgeEdgeCollisionsOld; };
  const vector<pair<VECTOR2,VECTOR2> >& edgeEdgeCoordinatesOld() const { return _edgeEdgeCoordinatesOld; };
  int totalCollisions() const { return _edgeEdgeCollisions.size(); };
  
  int totalVertices() const { return _vertices.size(); };
  virtual const int DOFs() const    { return _vertices.size() * 3 + _edges.size(); };
  const int totalEdges() const { return _totalEdges; };
  const REAL E() const       { return _E; };
  //const REAL nu() const      { return _nu; };
  const REAL G() const      { return _G; };
  const REAL density() const { return _density; };
  const REAL kStretching() const { return _E * M_PI * _radiusA * _radiusB; };
  bool& bendingForceFilterEnabled() { return _bendingForceFilterEnabled; };

  const REAL vertexMass(const int index)           { return _vertexMasses[index]; };

  // get the area of each edge
  const VECTOR voronoiAreas() const;
  
  // set the vertex displacements to these values exactly
  virtual void setDisplacement(const VECTOR& delta);
  virtual const VECTOR getDisplacement() const;
  const VECTOR getPositions() const;
  const VECTOR& thetas() const { return _thetas; };
  VECTOR& thetas() { return _thetas; };

  // recompute all the angle information, bending, twisting, frames, etc.
  void updateProperties();

  // bending energy
  virtual REAL computeBendingEnergy();
  virtual VECTOR computeBendingForces();
  virtual SPARSE_MATRIX computeBendingHessian();
  virtual SPARSE_MATRIX computeBendingClampedHessian();

  // just for verification purposes
  VECTOR computeBendingForcesFiniteDiff();
  MATRIX computeTwistingHessianFiniteDiff();

  // twisting energy
  virtual VECTOR computeTwistingForces() const;
  virtual SPARSE_MATRIX computeTwistingHessian() const;
  virtual SPARSE_MATRIX computeTwistingClampedHessian() const;

  // this is here just to debug twisting energies
  void cacheOldAngles();

  // compute stretching quantities
  virtual REAL computeStretchingEnergy(const STRAND::STRETCHING& hyperelastic) const;
  virtual VECTOR computeStretchingForces(const STRAND::STRETCHING& hyperelastic) const;
  virtual SPARSE_MATRIX computeStretchingHessian(const STRAND::STRETCHING& hyperelastic) const;
  virtual SPARSE_MATRIX computeStretchingClampedHessian(const STRAND::STRETCHING& hyperelastic) const;
  virtual SPARSE_MATRIX computeStretchingGaussNewtonHessian(const STRAND::STRETCHING& hyperelastic) const;

  // compute collision quantities
  virtual VECTOR computeEdgeEdgeCollisionForces(const ENERGY_12D& eeEnergy) const;
  virtual SPARSE_MATRIX computeEdgeEdgeCollisionClampedHessian(const ENERGY_12D& eeEnergy) const;

  // print out all the vertex positions and thetas
  void printState();

  // collision detection
  virtual void computeEdgeEdgeCollisions(const REAL& collisionEps, const bool verbose);

  // find the strand with the longest length
  const int longestStrand() const;

  // get the length of an individual strand
  const REAL strandLength(const unsigned int index) const;

  // output a single strand to an SOBJ file
  const bool writeStrand(const string& filename, const int index) const;

  // for scene replay visualization, allow collision pairs to be reset
  void setCollisionPairs(const vector<pair<int,int>>& eeColllsions);

  // File IO for strand-OBJ (SOBJ) files
  static bool readSOBJFile(const char* filename, vector<VECTOR3>& vertices,
                           vector<vector<int>>& indices);
  static bool writeSOBJFile(const char* filename,
                            const vector<VECTOR3>& vertices,
                            const vector<vector<int>>& indices);

protected:
  // generic initialization across multiple constructors
  virtual void initialize();

  // compute the reference twist currently induced by parallel transport only
  REAL currentReferenceTwist(const int i);

  void computeEdges();
  void computeReferenceDirectors();
  void computeTangents();
  void computeReferenceTwist();
  void computeCurvatureBinormals();
  void computeTanBinormals();
  void computeEdgeLengths();
  void computeVoronoiLengths();
  void computeMaterialDirectors();
  void computeVertexMasses();
  void computeSpaceParallel();

  void computeKappas();
  void computeTwists();
  void computeBs(const REAL& E, const REAL& radiusA, const REAL& radiusB);
  void computeRefVertexLengths();

  MATRIX11x2 computeGradKappa(const int bendIndex) const;
  pair<MATRIX11, MATRIX11> computeHessianKappa(const int bendIndex) const;
  VECTOR11 computeTanGradTwist(const int bendIndex) const; 
  MATRIX11 computeTanHessianTwist(const int bendIndex) const; 
  VECTOR11 computeDWAGradTwist(const int bendIndex) const; 

  // shared vector/matrix constructions so they don't have to be repeated
  // across clamped and unclamped Hessians
  SPARSE_MATRIX buildPerBendMatrix(const vector<MATRIX11>& perBendHessians) const;
  SPARSE_MATRIX buildPerEdgeMatrix(const vector<MATRIX6>& perEdgeHessians) const;
  virtual SPARSE_MATRIX buildEdgeEdgeMatrix(const vector<MATRIX12>& perEdgeHessians) const;
  VECTOR buildPerBendVector(const vector<VECTOR11>& perBendForces) const;
  VECTOR buildPerEdgeVector(const vector<VECTOR6>& perEdgeForces) const;
  VECTOR buildEdgeEdgeVector(const vector<VECTOR12>& perEdgeForces) const;

  // version that puts the edges at the end
  SPARSE_MATRIX buildPerBendMatrixEdgeEnd(const vector<MATRIX11>& perBendHessians) const;
  SPARSE_MATRIX buildPerEdgeMatrixEdgeEnd(const vector<MATRIX6>& perEdgeHessians) const;
  SPARSE_MATRIX buildEdgeEdgeMatrixEdgeEnd(const vector<MATRIX12>& perEdgeHessians) const;
  VECTOR buildPerBendVectorEdgeEnd(const vector<VECTOR11>& perBendForces) const;
  VECTOR buildPerEdgeVectorEdgeEnd(const vector<VECTOR6>& perEdgeForces) const;
  VECTOR buildEdgeEdgeVectorEdgeEnd(const vector<VECTOR12>& perEdgeForces) const;

  // version that interleaves the edge indices
  SPARSE_MATRIX buildPerBendMatrixInterleaved(const vector<MATRIX11>& perBendHessians) const;
  SPARSE_MATRIX buildPerEdgeMatrixInterleaved(const vector<MATRIX6>& perEdgeHessians) const;
  SPARSE_MATRIX buildEdgeEdgeMatrixInterleaved(const vector<MATRIX12>& perEdgeHessians) const;
  VECTOR buildPerBendVectorInterleaved(const vector<VECTOR11>& perBendForces) const;
  VECTOR buildPerEdgeVectorInterleaved(const vector<VECTOR6>& perEdgeForces) const;
  VECTOR buildEdgeEdgeVectorInterleaved(const vector<VECTOR12>& perEdgeForces) const;

  // version that puts the edges at the end
  void setDisplacementEdgeEnd(const VECTOR& delta);
  const VECTOR getDisplacementEdgeEnd() const;

  // version that interleaves the edges
  void setDisplacementInterleaved(const VECTOR& delta);
  const VECTOR getDisplacementInterleaved() const;

  // compute the vertex and edge indices inside the global vector
  void buildGlobalIndices();

  // Tangent-based forces
  VECTOR computeTanBendingForces();
  virtual SPARSE_MATRIX computeTanBendingHessian();
  virtual SPARSE_MATRIX computeTanBendingClampedHessian();
  virtual SPARSE_MATRIX computeBendingGaussNewtonHessian();
  virtual SPARSE_MATRIX computeTanBendingGaussNewtonHessian();
  
  // Tan-based twisting energy
  VECTOR computeTanTwistingForces() const;
  virtual SPARSE_MATRIX computeTanTwistingHessian() const;
  virtual SPARSE_MATRIX computeTanTwistingClampedHessian() const;
  virtual SPARSE_MATRIX computeTwistingGaussNewtonHessian() const;
  virtual SPARSE_MATRIX computeTanTwistingGaussNewtonHessian() const;
  
  // DreamWorks twisting energy
  VECTOR computeDWATwistingForces() const;
  virtual SPARSE_MATRIX computeDWATwistingHessian() const;
  virtual SPARSE_MATRIX computeDWATwistingClampedHessian() const;

  // are we putting the edge information at the end?
  bool _edgeEnd;

  // the core geometry
  vector<VECTOR3>    _vertices;
  vector<VECTOR3>    _restVertices;
  unsigned int       _totalVertices;
  unsigned int       _totalEdges;
  unsigned int       _totalBends;
  unsigned int       _totalStrands;
  unsigned int       _DOFs;         // 3 * _totalVertices + _totalEdges

  // old vertex positions, just for visualization
  vector<VECTOR3>    _verticesOld;
  VECTOR             _thetasOld;

  // vector of indices for each individual strand
  vector<vector<int> > _strandIndices;

  // deformation gradient
  VECTOR             _dmInvs;

  // indices of edges that experience stretch
  vector<VECTOR2I> _edgeIndices;  // _totalEdges
  vector<VECTOR3>  _edges;

  // indices of bends, i.e. pairs of edges, with a bend in between
  // the vertex ordering is:
  //                   //
  //         1         //
  //         o         //
  //        / \        //
  //       /   \       //
  //      /     \      //
  //   0 /       \ 2   //
  //    o         o    //
  //                   //
  vector<VECTOR3I> _bendVertices; // _totalBends
  // the edge ordering is:
  //                   //
  //                   //
  //         o         //
  //        / \        //
  //     0 /   \ 1     //
  //      /     \      //
  //     /       \     //
  //    o         o    //
  //                   //
  vector<VECTOR2I> _bendEdges;    // _totalBends

  // list of bools that say whether the vertex has an associated edge.
  // along one bend, the starting vertex owns the edge that comes right
  // after it, e.g. v0 owns e0, and v1 owns e1. 
  //                   //
  //        v1         //
  //         o         //
  //        / \        //
  //    e0 /   \ e1    //
  //      /     \      //
  //     /       \     //
  // v0 o         o    //
  //                   //
  VECTORI _vertexOwnsEdge;
  VECTORI _globalVertexIndices;
  VECTORI _globalEdgeIndices;
  VECTORI _globalStrandEnds;  // what index does each strand end on?

  // twist frames from previous timestep
  vector<VECTOR3> _tangentsOld;    // _totalEdges
  vector<VECTOR3> _directorOld1;   // _totalEdges
  vector<VECTOR3> _directorOld2;   // _totalEdges
  
  // rest curvatures
  vector<VECTOR2> _restKappas;      // _totalVertices

  vector<VECTOR3> _tangents;
  vector<VECTOR3> _director1s;      // _totalEdges
  vector<VECTOR3> _director2s;      // _totalEdges
  vector<VECTOR3> _material1s;      // _totalEdges
  vector<VECTOR3> _material2s;      // _totalEdges
  vector<MATRIX2> _perCornerKappas; // _totalVertices 

  // curvature binormals (kuravture?)
  vector<VECTOR3> _kbs;             // _totalBends
  vector<VECTOR2> _kappas;          // _totalBends
  vector<VECTOR2> _kappaBars;       // _totalBends
  vector<MATRIX2> _Bs;              // _totalBends
  VECTOR          _refVertexLengths;// _totalBends

  VECTOR          _referenceTwists;  // _totalBends
  VECTOR          _twists;           // _totalBends
  VECTOR          _undeformedTwists; // _totalBends
  VECTOR          _kts;              // _totalBends, twisting stiffness
  VECTOR          _thetas;           // _totalEdges
  VECTOR          _edgeLengths;      // _totalEdges
  VECTOR          _restEdgeLengths;  // _totalEdges
  VECTOR          _restBendLengths;  // _totalBends
  VECTOR          _restBendAngles;   // _totalBends
  VECTOR          _restTwistAngles;  // _totalEdges, though seems like it should be _totalVertices

  VECTOR          _voronoiLengths;    // _totalVertices
  VECTOR          _vertexMasses;      // _totalVertices

  // material parameters
  REAL _E;
  //REAL _nu; // let's prefer shear modulus to Poisson's ratio
  REAL _G;
  REAL _density;
  REAL _radiusA;
  REAL _radiusB;

  // list of edge-edge collision indices
  // first indexes into _surfaceEdges
  // second indexes into _surfaceEdges
  vector<pair<int, int> > _edgeEdgeCollisions;
  vector<pair<int, int> > _edgeEdgeCollisionsOld;

  // interpolation coordinates for edge-edge collisions
  vector<pair<VECTOR2, VECTOR2> > _edgeEdgeCoordinates;
  vector<pair<VECTOR2, VECTOR2> > _edgeEdgeCoordinatesOld;
  
  // are the edge-edge collisions still separate, or is there already a face-edge intersection?
  vector<bool> _edgeEdgeIntersections;
  
  // scaling term for edge-edge collision forces
  vector<REAL> _edgeEdgeCollisionAreas;
  
  // how close is considered to be in collision?
  //REAL _collisionEps;
  
  // mapping from edge index pairs to _edgeIndices
  map<pair<int, int>, int> _edgeHash;
  
  // which edge-edge collision force are we using?
  //VOLUME::EDGE_COLLISION* _edgeEdgeEnergy;

  // filter bending forces? Need to come up with a more forgiving energy here.
  bool _bendingForceFilterEnabled;
  REAL _bendingFilterThreshold;
};

}

#endif
