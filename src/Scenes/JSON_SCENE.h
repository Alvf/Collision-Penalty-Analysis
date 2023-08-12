#ifndef JSON_SCENE_H
#define JSON_SCENE_H

#include "SIMULATION_SCENE.h"

namespace HOBAK {

// Replaying simulation output from JSON file.
//
// Need this purely virtual class because volumes, shells, and strands
// all handle things slightly differently
class JSON_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override = 0;

  // don't do anything for building the scene - FILE_IO should be doing
  // the heavy lifting here
  virtual bool buildScene() override = 0;
  
  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override = 0;

  // jump to a specific frame
  virtual void jumpToFrame(const int frame) = 0;

  // in case _drawFrame changed, you can update the positions
  virtual void updatePositions() = 0;

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene() override = 0;

  // subtract off the translation and rotation from the main body
  // and draw that
  virtual void drawBodyCenteredScene(const VECTOR3& primitiveColor = VECTOR3(0.5, 0.5, 0.5), 
                                     const VECTOR3& outlineColor = VECTOR3(0,0,0)) = 0;
#endif

  const vector<VECTOR>& positions() const  { return _positions; };
  const vector<VECTOR>& velocities() const { return _velocities; };
  const vector<vector<pair<int,int> > >& vertexFaceCollisions() const { return _vertexFaceCollisions; };
  const vector<vector<pair<int,int> > >& edgeEdgeCollisions() const   { return _edgeEdgeCollisions; };
  vector<VECTOR>& velocities()             { return _velocities; };
  vector<VECTOR>& positions()              { return _positions; };
  vector<vector<pair<int,int> > >& vertexFaceCollisions() { return _vertexFaceCollisions; };
  vector<vector<pair<int,int> > >& edgeEdgeCollisions()   { return _edgeEdgeCollisions; };
  int& drawFrame()                         { return _drawFrame; };
  const int totalFrames() const            { return _positions.size(); };

  // have a setter here for data hiding
  void setInitialA(const MATRIX3& A)           { _initialA = A; };
  void setInitialTranslation(const VECTOR3& t) { _initialTranslation = t; };

protected:
  vector<VECTOR> _positions;
  vector<VECTOR> _velocities;
  vector<vector<pair<int,int> > > _vertexFaceCollisions;
  vector<vector<pair<int,int> > > _edgeEdgeCollisions;
  int _drawFrame;
};

}

#endif
