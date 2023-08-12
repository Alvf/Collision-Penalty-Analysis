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
#ifndef TET_MESH_H
#define TET_MESH_H

#include "SETTINGS.h"
#include "util/BLOCK_DIAGONAL_MATRIX3.h"
#include "util/BLOCK_SPARSE_MATRIX3.h"
#include "Hyperelastic/Volume/HYPERELASTIC.h"
#include "Collision/ENERGY_12D.h"
// included so we can collide against a strand mesh
#include "STRAND_MESH.h"
#include "Damping/Volume/DAMPING.h"
#include "util/BLOCK_SPARSE_MATRIX3.h"

#include <map>
#include <vector>

namespace HOBAK {

using namespace std;

class TET_MESH
{
public:
  TET_MESH() = default;
  TET_MESH(const vector<VECTOR3>& restVertices, 
           const vector<VECTOR4I>& tets);
  virtual ~TET_MESH();

  // accessors
  const vector<VECTOR3>& vertices() const     { return _vertices; };
  vector<VECTOR3>& vertices()                 { return _vertices; };
  const vector<VECTOR3>& restVertices() const { return _restVertices; };
  vector<VECTOR3>& restVertices()             { return _restVertices; };
  const vector<VECTOR4I>& tets() const        { return _tets; };
  vector<VECTOR4I>& tets()                    { return _tets; };
  const vector<VECTOR4I>& vertexFaceCollisionTets() const    { return _vertexFaceCollisionTets; };
  vector<VECTOR4I>& vertexFaceCollisionTets()                { return _vertexFaceCollisionTets; };
  const vector<REAL>& restOneRingVolumes() const   { return _restOneRingVolumes; };
  const vector<REAL>& restOneRingAreas() const   { return _restOneRingAreas; };
  const VECTOR& restEdgeAreas() const   { return _restEdgeAreas; };
  vector<REAL>& restOneRingVolumes()               { return _restOneRingVolumes; };
  const VECTOR3& vertex(const int index) const     { return _vertices[index]; };
  VECTOR3& vertex(const int index)                 { return _vertices[index]; };
  const VECTOR3& restVertex(const int index) const { return _restVertices[index]; };
  const VECTOR4I& tet(const int tetIndex) const    { return _tets[tetIndex]; };
  const vector<int>& surfaceTets() const           { return _surfaceTets; };
  const vector<int>& surfaceVertices() const       { return _surfaceVertices; };
  const vector<int>& toSurfaceVertices() const       { return _toSurfaceVertices; };
  const vector<VECTOR3I>& surfaceTriangles() const { return _surfaceTriangles; };
  const vector<VECTOR3I>& surfaceTrianglesIntoSurfaceVertices() const { return _surfaceTrianglesIntoSurfaceVertices; };
  const vector<VECTOR2I>& surfaceEdges() const     { return _surfaceEdges; };
  const vector<pair<int,int> >& vertexFaceCollisions() const { return _vertexFaceCollisions; };
  const vector<pair<int,int> >& edgeEdgeCollisions() const   { return _edgeEdgeCollisions; };
  const vector<bool>& edgeEdgeIntersections() const          { return _edgeEdgeIntersections; };
  const vector<pair<VECTOR2,VECTOR2> >& edgeEdgeCoordinates() const   { return _edgeEdgeCoordinates; };
  const vector<REAL>& surfaceTriangleAreas() const { return _surfaceTriangleAreas; };
  const vector<VECTOR3I>& surfaceTriangleNeighbors() const { return _surfaceTriangleNeighbors; };
  const vector<bool>& invertedVertices() const {return _invertedVertices; };
  map<int, int>& volumeToSurfaceID() {return _volumeToSurfaceID; };
  const vector<VECTOR2I> surfaceEdgeTriangleNeighbors() const  {return _surfaceEdgeTriangleNeighbors; };
  const vector<bool> visibility() const {return _visibility; };
  vector<bool>& visibility() {return _visibility; };
  const vector<VECTOR3I> visibleSurfaceTriangles() const {return _visibleSurfaceTriangles; };
  vector<VECTOR3I>& visibleSurfaceTriangles()            {return _visibleSurfaceTriangles; };
  const vector<VECTOR3I> visibleSurfaceTrianglesIntoSurfaceVertices() const {return _visibleSurfaceTrianglesIntoSurfaceVertices; };
  vector<VECTOR3I>& visibleSurfaceTrianglesIntoSurfaceVertices()            {return _visibleSurfaceTrianglesIntoSurfaceVertices; };

  int totalVertices() const { return _vertices.size(); };
  const int DOFs() const    { return _vertices.size() * 3; };

  // get deformation gradient, and its SVD
  MATRIX3 computeF(const int tetIndex) const;
  void computeFs();
  void computeFdots(const VECTOR& velocity);
  void computeSVDs();

  // get volume-weighted global translation
  VECTOR3 getTranslation() const;
  
  // get volume-weighted global translation, for the rest state
  VECTOR3 getRestTranslation() const;
  
  // get Procrustes-style global rotation
  MATRIX3 getRotation() const;

  // get the current displacement in vector form
  VECTOR getDisplacement() const;
  
  // set the vertex displacements to these values exactly
  void setDisplacement(const VECTOR& delta);

  // set the vertex positions directly exactly
  void setPositions(const VECTOR& positions);

  // add the following deltas to the positions
  void addDisplacement(const VECTOR& delta);

  // set collision pairs (for replays)
  void setCollisionPairs(const vector<pair<int,int> >& vertexFace, 
                         const vector<pair<int,int> >& edgeEdge);

  // compute hyperelastic quantities
  REAL computeHyperelasticEnergy(const VOLUME::HYPERELASTIC& hyperelastic) const;
  VECTOR computeHyperelasticForces(const VOLUME::HYPERELASTIC& hyperelastic) const;
  virtual SPARSE_MATRIX computeHyperelasticClampedHessian(const VOLUME::HYPERELASTIC& hyperelastic) const;
  virtual SPARSE_MATRIX computeHyperelasticHessian(const VOLUME::HYPERELASTIC& hyperelastic) const;
  virtual BLOCK_SPARSE_MATRIX3 computeBlockHyperelasticClampedHessian(const VOLUME::HYPERELASTIC& hyperelastic) const;

  // compute damping quantities
  VECTOR computeDampingForces(const VOLUME::DAMPING& damping) const;
  virtual SPARSE_MATRIX computeDampingHessian(const VOLUME::DAMPING& damping) const;

  // compute collision quantities
  //VECTOR computeCollisionForces(const REAL& collisionStiffness) const;
  //SPARSE_MATRIX computeCollisionClampedHessian(const REAL& collisionStiffness) const;

  // compute x-based collision quantities
  VECTOR computeVertexFaceCollisionForces(const ENERGY_12D& vertexFaceEnergy) const;
  SPARSE_MATRIX computeVertexFaceCollisionClampedHessian(const ENERGY_12D& vertexFaceEnergy) const;
  //VECTOR computeEberleDampingForces(const VECTOR& velocity, const REAL& integratorConstant,
  //                                  const REAL& collisionStiffness) const;
  //SPARSE_MATRIX computeEberleClampedHessian(const VECTOR& velocity, const REAL& integratorConstant,
  //                                          const REAL& collisionStiffness) const;
  REAL computeEdgeEdgeCollisionEnergy(const ENERGY_12D& edgeEdgeEnergy) const;
  VECTOR computeEdgeEdgeCollisionForces(const ENERGY_12D& edgeEdgeEnergy) const;
  //VECTOR computeEdgeEdgeDampingForces(const VECTOR& velocity, const REAL& integratorConstant,
  //                                    const REAL& collisionStiffness) const;
  SPARSE_MATRIX computeEdgeEdgeCollisionClampedHessian(const ENERGY_12D& edgeEdgeEnergy) const;
  //SPARSE_MATRIX computeEdgeEdgeDampingClampedHessian(const VECTOR& velocity, const REAL& integratorConstant,
  //                                                   const REAL& collisionStiffness) const;

  // compuate elastic and damping forces at the same time
  virtual VECTOR computeInternalForces(const VOLUME::HYPERELASTIC& hyperelastic,
                                       const VOLUME::DAMPING& damping) const;

  // get the bounding box for the current mesh
  void getBoundingBox(VECTOR3& mins, VECTOR3& maxs) const;

  // find all the vertex-face collision pairs, using the InFaceRegion test
  virtual void computeVertexFaceCollisions(const REAL& collisionEps);

  // find all the edge-edge collision pairs
  virtual void computeEdgeEdgeCollisions(const REAL& collisionEps);

  // collisions with strands
  virtual void computeVertexFaceCollisionsWithStrands(const STRAND_MESH& strandMesh,
                                                      const ENERGY_12D& vertexFaceEnergy, 
                                                      vector<pair<int, int> >& collisions);
  virtual void computeEdgeEdgeCollisionsWithStrands(const STRAND_MESH& strandMesh, 
                                                    const ENERGY_12D& edgeEdgeEnergy, 
                                                    vector<pair<int, int> >& collisions,
                                                    vector<bool>& intersections,
                                                    vector<pair<VECTOR2, VECTOR2> >& coordinates,
                                                    vector<REAL>& areas) const;

  // debug edge-edge collisions, load up some specific pairs
  //void computeEdgeEdgeCollisionsDebug();

  // based on vertex-face collision pairs, build "collision tets"
  void buildVertexFaceCollisionTets(const VECTOR& velocity);

  // based on vertex-face collision pairs, build "collision tets" for strand volume collisions
  void buildVertexFaceCollisionTetsWithStrands(const STRAND_MESH& strandMesh, 
                                               const vector<pair<int, int> >& collisions,
                                               vector<VECTOR4I>& tets,
                                               vector<REAL>& areas) const;
  
  // write out the surface to OBJ triangle mesh
  static bool writeSurfaceToObj(const string& filename, const TET_MESH& tetMesh);

  // read in OBJ-style tet mesh file
  static bool readTobjFile(const string& filename, 
                           vector<VECTOR3>& vertices, vector<VECTOR4I>& tets,
                           bool addingOn = false, const VECTOR3& translate = VECTOR3(0,0,0));

  // write out OBJ-style tet mesh file, the bool tells whether to write the
  // rest or deformed vertices
  static bool writeTobjFile(const string& filename, const TET_MESH& tetMesh, const bool restVertices);

  // normalize vertices so that they're in a unit box, centered at (0.5, 0.5, 0.5)
  static vector<VECTOR3> normalizeVertices(const vector<VECTOR3>& vertices);

  // compute distance between a point and triangle
  static REAL pointTriangleDistance(const VECTOR3& v0, const VECTOR3& v1, 
                                    const VECTOR3& v2, const VECTOR3& v);

  // see if the projection of v onto the plane of v0,v1,v2 is inside the triangle
  // formed by v0,v1,v2
  static bool pointProjectsInsideTriangle(const VECTOR3& v0, const VECTOR3& v1, 
                                          const VECTOR3& v2, const VECTOR3& v);

  // compute the dihedral angle between two surface faces
  REAL surfaceFaceDihedralAngle(const int surfaceID0, const int surfaceID1) const;

  // see if a current surface triangle has been crushed to degeneracy
  bool surfaceTriangleIsDegenerate(const int surfaceTriangleID);

  // compute a triangle area
  static REAL triangleArea(const vector<VECTOR3>& triangle);

  // see if the vertex is inside the collision cell described in 
  // Chapter 11: Collision Processing, [Kim and Eberle 2020]
  bool insideCollisionCell(const int surfaceTriangleID, const VECTOR3& vertex);
protected:
  // compute the volume of a tet
  static REAL computeTetVolume(const vector<VECTOR3>& tetVertices);

  // compute volumes for tets -- works for rest and deformed, just
  // pass it _restVertices or _vertices
  vector<REAL> computeTetVolumes(const vector<VECTOR3>& vertices);

  // compute volumes for each vertex one-ring -- works for rest
  // and deformed, just pass it _restVertices or _vertices
  vector<REAL> computeOneRingVolumes(const vector<VECTOR3>& vertices);

  // compute material inverses for deformation gradient
  vector<MATRIX3> computeDmInvs();

  // compute change-of-basis from deformation gradient F to positions, x
  vector<MATRIX9x12> computePFpxs();

  // find what's on the surface
  void computeSurfaceVertices();
  void computeSurfaceTriangles();
  void computeSurfaceEdges();
  void computeSurfaceAreas();
  void computeSurfaceTriangleNeighbors();
  void computeSurfaceEdgeTriangleNeighbors();

  // (DEACTIVATED) find out how close all the edges are initially
  void computeEdgeEdgeRestDistance();

  // compute a triangle area
  // static REAL triangleArea(const vector<VECTOR3>& triangle);

  // get the normal to a plane, specified by three points
  static VECTOR3 planeNormal(const vector<VECTOR3>& plane);

  // project point onto plane, specific by three points
  static VECTOR3 pointPlaneProjection(const vector<VECTOR3>& plane, const VECTOR3& point);

  // see if the vertex is inside the collision cell described in 
  // Chapter 11: Collision Processing, [Kim and Eberle 2020]
  // bool insideCollisionCell(const int surfaceTriangleID, const VECTOR3& vertex);

  // compute distance to collision cell wall, where positive means inside
  // and negative means outside
  REAL distanceToCollisionCellWall(const int surfaceTriangleID, const VECTOR3& vertex);

  // compute whether one vertex is inside the vertex one right of another
  void computeSurfaceVertexOneRings();

  // are these two surface triangles neighbors?
  bool areSurfaceTriangleNeighbors(const int id0, const int id1) const;

  // build a consistent tet/flap ordering from two surface triangles
  VECTOR4I buildSurfaceFlap(const int surfaceID0, const int surfaceID1) const;

  // compute the normal of the surface triangle at _surfaceTriangles[triangleID];
  VECTOR3 surfaceTriangleNormal(const int triangleID) const;

  // see if a current surface triangle has been crushed to degeneracy
  // bool surfaceTriangleIsDegenerate(const int surfaceTriangleID);

  // compute which vertices are attached to inverted tets
  void computeInvertedVertices();

  // the core geometry
  vector<VECTOR3>    _vertices;
  vector<VECTOR3>    _restVertices;
  vector<VECTOR4I>   _tets;

  // indexes into surfaceTriangles
  vector<bool>       _visibility;
  vector<VECTOR3I>   _visibleSurfaceTriangles;
  vector<VECTOR3I>   _visibleSurfaceTrianglesIntoSurfaceVertices;

  // volumes, computed by computeTetVolumes and computeOneRingVolumes
  //
  // TODO: why aren't these just VECTORs?
  vector<REAL>  _restTetVolumes;
  vector<REAL>  _restOneRingVolumes;
  vector<REAL>  _restOneRingAreas;
  VECTOR        _restEdgeAreas;

  // support for computing deformation gradient F
  vector<MATRIX3>    _DmInvs;

  // change-of-basis to go from deformation gradient (F) to positions (x)
  vector<MATRIX9x12> _pFpxs;

  // deformation gradients, and their SVDs
  vector<MATRIX3> _Fs;
  vector<MATRIX3> _Us;
  vector<VECTOR3> _Sigmas;
  vector<MATRIX3> _Vs;

  // velocity gradients
  vector<MATRIX3> _Fdots;

  // list of tets that are on the surface
  vector<int> _surfaceTets;

  // list of triangles that are on the surface
  // each triplet is ordered counter-clockwise, facing outwards
  // the VECTOR3I indexes into _vertices
  vector<VECTOR3I> _surfaceTriangles;
  // the VECTOR3I indexes into _surfaceVertices
  vector<VECTOR3I> _surfaceTrianglesIntoSurfaceVertices;
  vector<REAL> _surfaceTriangleAreas;

  // for each surface triangle, what's the index of the neighboring triangles?
  vector<VECTOR3I> _surfaceTriangleNeighbors;

  // list of edges on the surface
  // each pair is in sorted order, and index into _vertices
  vector<VECTOR2I> _surfaceEdges;

  // list of vertices that are on the surface
  // indexes into _vertices
  vector<int> _surfaceVertices;

  // map the vertex index to surface vertices index
  vector<int> _toSurfaceVertices;

  // for each _surfaceEdges, what are the one or two neighboring triangles
  // in _surfaceTriangles?
  vector<VECTOR2I> _surfaceEdgeTriangleNeighbors;

  // for each pair of _surfaceEdges, what collisionEps should we use? If they started
  // out closer than collisionEps, then we need to set a smaller tolerance.
  //
  // the entries in the pair<int,int> are:
  //   unsigned int flat = edge[0] + edge[1] * _surfaceEdges.size();
  map<pair<unsigned int, unsigned int>, REAL> _edgeEdgeRestDistance;

  // how close is considered to be in collision?
  //REAL _collisionEps;

  // list of vertex-face collisions
  // first indexes into _vertices
  // second indexes into _surfaceTriangles
  vector<pair<int, int> > _vertexFaceCollisions;

  // list of edge-edge collision indices
  // first indexes into _surfaceEdges
  // second indexes into _surfaceEdges
  vector<pair<int, int> > _edgeEdgeCollisions;

  // interpolation coordinates for edge-edge collisions
  vector<pair<VECTOR2, VECTOR2> > _edgeEdgeCoordinates;
  
  // are the edge-edge collisions still separate, or is there already a face-edge intersection?
  vector<bool> _edgeEdgeIntersections;

  // list of collisionEps for each _edgeEdgeIntersections. Usually this will be the default,
  // but if it started out closer than that, we use that instead
  //vector<REAL> _edgeEdgeCollisionEps;

  // list of "collision tets" formed by vertex-face pairs
  vector<VECTOR4I> _vertexFaceCollisionTets;

  // list of "collision tets" formed by edge-edge pairs
  //vector<VECTOR4I> _edgeEdgeCollisionsTets;

  // DEBUG: see if the collision tet exists already
  //map<pair<int, int>, int> _vertexFaceCollisionTetsHash;
  
  // scaling term for vertex-face collision forces
  vector<REAL> _vertexFaceCollisionAreas;

  // scaling term for edge-edge collision forces
  vector<REAL> _edgeEdgeCollisionAreas;

  // convert tet mesh vertexID into a surface mesh vertexID
  // convert into into _vertices into index into _surfaceVertices
  map<int, int> _volumeToSurfaceID;

  // have your computed the SVDs since the last time you computed F?
  bool _svdsComputed;

  // see if two indices in _vertices (in sorted order)
  // are within the one ring of each other
  map<pair<int,int>, bool> _insideSurfaceVertexOneRing;

  // which vertices are inverted?
  vector<bool> _invertedVertices;
};

}

#endif
