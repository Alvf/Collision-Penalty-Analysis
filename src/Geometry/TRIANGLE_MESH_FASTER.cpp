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
#include "TRIANGLE_MESH_FASTER.h"
#include "MATRIX_UTIL.h"
#include <iostream>
#include <float.h>
#include "util/TIMER.h"
#include "Collision/COLLISION_UTIL.h"
#include "LINE_INTERSECT.h"

// #define VERY_VERBOSE 1

namespace HOBAK {
using namespace std;
TRIANGLE_MESH_FASTER::TRIANGLE_MESH_FASTER(const vector<VECTOR3>& restVertices, const vector<VECTOR3I>& triangles) : 
    TRIANGLE_MESH(restVertices, triangles),
    _aabbTreeTriangles(_vertices, &_triangles),
    _aabbTreeEdges(_vertices, &_edges)
{
  // mapping from edge index pairs to _surfaceEdges
  for (unsigned int x = 0; x < _edges.size(); x++)
  {
    pair<int,int> edge(_edges[x][0], _edges[x][1]);
    _edgeHash[edge] = x;
  }
}

// Finding edge-edge collisions with AABB tree helping the inner loop
void TRIANGLE_MESH_FASTER::computeEdgeEdgeCollisions()
{
  TIMER functionTimer(string("TRIANGLE_MESH_FASTER::") + string(__FUNCTION__));
  _edgeEdgeCollisions.clear();
  _edgeEdgeIntersections.clear();
  //_edgeEdgeCollisionEps.clear();
  _edgeEdgeCoordinates.clear();
  _edgeEdgeCollisionAreas.clear();
  
  _aabbTreeEdges.refit();

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
    vector<int> coarseEdges;
    _aabbTreeEdges.nearbyEdges(_edges[x],_collisionEps,coarseEdges);
    
    // find the closest other edge
    for (unsigned int y = 0; y < coarseEdges.size(); y++)
    {

      // skip if index is smaller -- don't want to double count nearby edges
      // (a,b) and (b,a)
      if ((unsigned int)coarseEdges[y] < x) continue;

      const VECTOR2I innerEdge = _edges[coarseEdges[y]];
      // if they share a vertex, skip it
      if ((outerEdge[0] == innerEdge[0]) || (outerEdge[0] == innerEdge[1]) ||
          (outerEdge[1] == innerEdge[0]) || (outerEdge[1] == innerEdge[1]))
        continue;

      const VECTOR3& v2 = _vertices[innerEdge[0]];
      const VECTOR3& v3 = _vertices[innerEdge[1]];

      /*
      // call the geometric test
      VECTOR3 innerPoint, outerPoint, midpoint, normal;
      bool intersect;
      IntersectLineSegments(v0[0], v0[1], v0[2], v1[0], v1[1], v1[2],
                            v2[0], v2[1], v2[2], v3[0], v3[1], v3[2],
                            false, 1e-8,
                            outerPoint[0], outerPoint[1], outerPoint[2],
                            innerPoint[0], innerPoint[1], innerPoint[2],
                            midpoint[0], midpoint[1], midpoint[2],
                            normal[0], normal[1], normal[2], intersect);
                            */
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
      closestEdge = coarseEdges[y];

      aClosest = a;
      bClosest = b;
    }

    // if nothing was close, move on
    if (closestEdge == -1) continue;

    /*
    // retrieve the eps of the closest edge
    const VECTOR2I innerEdge = _edges[closestEdge];
    const unsigned int innerFlat = innerEdge[0] + innerEdge[1] * _edges.size();
    const pair<unsigned int, unsigned int> edgeEdge(innerFlat, outerFlat);

    // it exists, right?
    assert(_edgeEdgeRestDistance.find(edgeEdge) != _edgeEdgeRestDistance.end());
    //const REAL collisionEps = _edgeEdgeRestDistance[edgeEdge];
    const REAL eeDistance = _edgeEdgeRestDistance[edgeEdge];
    const REAL collisionEps = (eeDistance < _collisionEps) ? eeDistance : _collisionEps;
    */

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

    /*
    const VECTOR2I innerEdge = _edges[closestEdge];
    pair<int, int> lookup;
    lookup.first = outerEdge[0];
    lookup.second = innerEdge[0];
    if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end()) continue;
    lookup.first = outerEdge[0];
    lookup.second = innerEdge[1];
    if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end()) continue;
    lookup.first = outerEdge[1];
    lookup.second = innerEdge[0];
    if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end()) continue;
    lookup.first = outerEdge[1];
    lookup.second = innerEdge[1];
    if (_insideSurfaceVertexOneRing.find(lookup) != _insideSurfaceVertexOneRing.end()) continue;
    */

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
      const REAL xArea = _restEdgeAreas[_edgeHash[outerPair]];
      const REAL closestArea = _restEdgeAreas[_edgeHash[innerPair]];
      _edgeEdgeCollisionAreas.push_back(xArea + closestArea);

      // find out if they are penetrating
      vector<VECTOR3> edge(2);
      edge[0] = v0;
      edge[1] = v1;

      // get the adjacent triangles of the *other* edge
      VECTOR2I adjacentTriangles = _edgeTriangleNeighbors[_edgeHash[innerPair]];

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
      //_edgeEdgeCollisionEps.push_back(_collisionEps);

      // TODO: for completeness, should probably test the other edges against the other
      // pair, just in case we're looking at a degenerate case. In general, seems redundant.
    }
  }
  assert(_edgeEdgeCollisions.size() == _edgeEdgeCoordinates.size());

#if VERY_VERBOSE
  if (_edgeEdgeCollisions.size() > 0)
    cout << " Collision area array size: " << _edgeEdgeCollisionAreas.size() << endl;
    cout << " Intersections array size: " << _edgeEdgeIntersections.size() << endl;
    cout << " Found " << _edgeEdgeCollisions.size() << " edge-edge collisions " << endl;

  cout << " pairs: " << endl;
  for (unsigned int x = 0; x < _edgeEdgeCollisions.size(); x++)
  {
    const pair<int,int> collision = _edgeEdgeCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")" << endl;
  }
#endif
}

// Finding vertex-face collisions with AABB tree helping the inner loop
void TRIANGLE_MESH_FASTER::computeVertexFaceCollisions()
{
  TIMER functionTimer(string("TRIANGLE_MESH_FASTER::") + string(__FUNCTION__));
  // if a vertex is part of an inverted triangle, don't have it participate 
  // in a self-collision. That triangle needs to get its house in order 
  // before it starts bossing around a surface face. Not checking for 
  // this causes faces to get horribly tangled in inverted configurations.
  computeInvertedVertices();

  _vertexFaceCollisions.clear();
  _aabbTreeTriangles.refit();

  for (unsigned int x = 0; x < _vertices.size(); x++)
  {
    // if the vertex is involved in an inverted triangle, give up
    if (_invertedVertices[x]) 
      continue;

    const VECTOR3& surfaceVertex = _vertices[x];
    vector<int> broadPhaseFaces;

    _aabbTreeTriangles.nearbyTriangles(surfaceVertex, _collisionEps, broadPhaseFaces);

    // find the close triangles
    for (unsigned int y = 0; y < broadPhaseFaces.size(); y++)
    {
      int faceID = broadPhaseFaces[y];
      
      // if the surface triangle is so small the normal could be degenerate, skip it
      if (surfaceTriangleIsDegenerate(faceID))
        continue;

      const VECTOR3I& t = _triangles[faceID];

      // if it's an inverted face, move on
      if (_invertedVertices[t[0]] && _invertedVertices[t[1]] && _invertedVertices[t[2]])
        continue;

      // if this triangle is in the one-ring of the current vertex, skip it
      if (t[0] == x || t[1] == x || t[2] == x) continue;

      const REAL distance = pointTriangleDistance(_vertices[t[0]], _vertices[t[1]],
                                                  _vertices[t[2]], surfaceVertex);

      if (distance < _collisionEps)
      {
        // if the point, projected onto the face's plane, is inside the face,
        // then record the collision now
        if (pointProjectsInsideTriangle(_vertices[t[0]], _vertices[t[1]], 
                                        _vertices[t[2]], surfaceVertex))
        {
          pair<int,int> collision(x, faceID);
          _vertexFaceCollisions.push_back(collision);
          continue;
        }
        if (insideCollisionCell(faceID, surfaceVertex))
        {
          pair<int,int> collision(x, faceID);
          _vertexFaceCollisions.push_back(collision);
        }
      }
    }
  }
#if VERY_VERBOSE
  if (_vertexFaceCollisions.size() > 0)
    cout << " Found " << _vertexFaceCollisions.size() << " vertex-face collisions " << endl;

  cout << " pairs: " << endl;
  for (unsigned int x = 0; x < _vertexFaceCollisions.size(); x++)
  {
    const pair<int,int> collision = _vertexFaceCollisions[x];
    cout << "(" << collision.first << ", " << collision.second << ")" << endl;
  }
#endif 
}

} // HOBAK
