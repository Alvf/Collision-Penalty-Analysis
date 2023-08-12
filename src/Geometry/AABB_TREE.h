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
#ifndef AABB_TREE_H
#define AABB_TREE_H

#include "KINEMATIC_SHAPE.h"

namespace HOBAK {

using namespace std;

// tree nodes for AABB_TREE
struct AABB_NODE
{
  AABB_NODE(const vector<int>& inputPrimitiveIndices, 
            const VECTOR3& inputMins, const VECTOR3& inputMaxs, 
            const int& inputDepth) : 
    mins(inputMins), maxs(inputMaxs), 
    primitiveIndices(inputPrimitiveIndices), depth(inputDepth)
  { };

  // what are the axis-aligned bounds of the box?
  VECTOR3 mins, maxs;

  // if we're an interior node, here are the children
  AABB_NODE* child[2] = {NULL, NULL};

  // indices of triangles in _surfaceTriangles, or edges in _surfaceEdges
  // enclosed by this AABB
  //
  // this list is only populated for leaf nodes
  vector<int> primitiveIndices;

  // track the depth we're at
  int depth;
};

/////////////////////////////////////////////////////////////////////////////////////////////
// Collision detection structure for vertex-triangle and edge-edge collision detection
//
// Nothing fancy here, just a search tree that gets the job done. There are online libarries
// that do something similar, but nothing that did quite what I want. Right now collision
// detection is consuming 90% of the running time, so anything is better than nothing.
/////////////////////////////////////////////////////////////////////////////////////////////
class AABB_TREE
{
public:
  AABB_TREE(const vector<VECTOR3>& vertices, const vector<VECTOR3I>* surfaceTriangles);
  AABB_TREE(const vector<VECTOR3>& vertices, const vector<VECTOR2I>* surfaceEdges);
	~AABB_TREE();

  // return a list of potential triangles nearby a vertex, subject to a distance threshold
  void nearbyTriangles(const VECTOR3& vertex, const REAL& eps, vector<int>& faces) const;
  
  // return a list of potential triangles nearby an edge, subject to a distance threshold
  void nearbyTriangles(const VECTOR2I& edge, const REAL& eps, vector<int>& faces) const;

  // return a list of potential nearby edges, subject to a distance threshold
  void nearbyEdges(const VECTOR2I& edge, const REAL& eps, vector<int>& faces) const;

  // return a list of potential nearby edges, subject to a distance threshold
  // here for debugging purposes, similar to the edge-based one, but doesn't
  // need to reference points in _vertices
  void nearbyTriangles(const VECTOR3& mins, const VECTOR3& maxs, const REAL& eps, vector<int>& faces) const;

  // get the root node
  const AABB_NODE& root() const { return *_root; };

  // refit the bounding boxes, presumably because the vertices moved
  void refit();

private:
  // build the tree for triangles!
  void buildTriangleRoot();
  
  // build the tree for edges!
  void buildEdgeRoot();

  // let's clean up after ourselves!
  void deleteTree(AABB_NODE* node);

  // given a triangle node, let's build it's children
  void buildTriangleChildren(AABB_NODE* node, const int depth);
  
  // given an edge node, let's build it's children
  void buildEdgeChildren(AABB_NODE* node, const int depth);

  // cut the a list of triangles into two child lists
  void buildTriangleChildLists(const REAL& cuttingPlane, const int& axis,
                               const vector<int>& triangleIndices,
                               vector<int>& childList0, vector<int>& childList1);

  // cut the a list of triangles into two child lists
  void buildEdgeChildLists(const REAL& cuttingPlane, const int& axis,
                           const vector<int>& edgeIndices,
                           vector<int>& childList0, vector<int>& childList1);

  // given a list of indices into _surfaceTriangles, find the 
  // bounding box
  void findTriangleBoundingBox(const vector<int>& triangleIndices, VECTOR3& mins, VECTOR3& maxs) const;

  // given a list of indices into _surfaceEdges, find the 
  // bounding box
  void findEdgeBoundingBox(const vector<int>& edgeIndices, VECTOR3& mins, VECTOR3& maxs) const;

  // recursively refit triangles
  void refitTriangles(AABB_NODE* node);

  // recursively refit edges
  void refitEdges(AABB_NODE* node);

  // return a list of potential triangles nearby a vertex, subject to a distance threshold
  void nearbyEdges(const AABB_NODE* node, const VECTOR2I& edge, 
                   const REAL& eps, vector<int>& edges) const;

  // return a list of potential triangles nearby a vertex, subject to a distance threshold
  void nearbyTriangles(const AABB_NODE* node, const VECTOR3& vertex, 
                       const REAL& eps, vector<int>& faces) const;
  
  // return a list of potential triangles nearby a box specified by min and max, 
  // subject to a distance threshold
  void nearbyTriangles(const AABB_NODE* node, const VECTOR3& mins, const VECTOR3& maxs, 
                       const REAL& eps, vector<int>& faces) const;

  // are we inside this AABB, subject to the distance eps?
  bool insideAABB(const AABB_NODE* node, const VECTOR3& vertex, const REAL& eps) const;
 
  // are these two AABBs overlapping, subject to the distance eps?
  bool overlappingAABBs(const AABB_NODE* node, const VECTOR3& mins, 
                        const VECTOR3& maxs, const REAL& eps) const;

  // do the AABB of this node and the AABB of the edge overlap?
  bool overlappingAABBs(const AABB_NODE* node, const VECTOR2I& edge, const REAL& eps) const;

  const vector<VECTOR3>& _vertices;

  // const pointer to the surface triangles in the tet mesh
  // make this a pointer so it can be NULL, in case we're building an edge tree
  const vector<VECTOR3I>* _surfaceTriangles;

  // const pointer to the surface edges in the tet mesh
  // make this a pointer so it can be NULL, in case we're building an edge tree
  const vector<VECTOR2I>* _surfaceEdges;

  AABB_NODE* _root;
};

}

#endif
