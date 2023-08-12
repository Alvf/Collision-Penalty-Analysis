#ifndef JSON_STRAND_SCENE_H
#define JSON_STRAND_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"

namespace HOBAK {

// Replaying simulation output from JSON file.
class JSON_STRAND_SCENE : public JSON_SCENE
{
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Replaying JSON output from strand simulation" << endl;
    cout << "=====================================================================" << endl;
  }

  // don't do anything for building the scene - FILE_IO should be doing
  // the heavy lifting here
  virtual bool buildScene() override
  {
    vector<VECTOR3> vertices;
    vector<vector<int>> indices;
    bool success = STRAND_MESH::readSOBJFile(_strandMeshFilename.c_str(), vertices, indices);

    if (!success)
    {
      cout << " Failed to open file " << _strandMeshFilename.c_str() << endl;
      return false;
    }

    // TODO: is this support needed?
    //if (_normalizedVertices)
    //  vertices = TRIANGLE_MESH::normalizeVertices(vertices);

    // apply any initial transforms the scene requested
    for (unsigned int x = 0; x < vertices.size(); x++)
      vertices[x] = _initialA * vertices[x] + _initialTranslation;

    _strandMesh = new TET_WISP_MESH(vertices, indices,
                                    _E, _G, _density, _baseRadius, _tipRadius);

    _drawFrame = 0;
    _drawFeature = true;
    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) 
  {
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    
    _strandMesh->setDisplacement(_positions[_drawFrame]);
    _strandMesh->setCollisionPairs(_edgeEdgeCollisions[_drawFrame]);

    if (verbose)
      cout << " Setting to frame " << _drawFrame << endl;
    _drawFrame++;
  };

  // jump to a specific frame
  void jumpToFrame(const int frame) 
  {
    _drawFrame = frame;
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    
    _strandMesh->setDisplacement(_positions[_drawFrame]);
    _strandMesh->setCollisionPairs(_edgeEdgeCollisions[_drawFrame]);

    cout << " Jumping to frame " << _drawFrame << endl;
  };

  // in case _drawFrame changed, you can update the positions
  void updatePositions()
  {
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    _strandMesh->setDisplacement(_positions[_drawFrame]);
    _strandMesh->setCollisionPairs(_edgeEdgeCollisions[_drawFrame]);
  }

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene() override
  {
    glEnable(GL_DEPTH_TEST);

    drawWispMesh(*_strandMesh, -1, false);

    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
  };

  // subtract off the translation and rotation from the main body
  // and draw that
  void drawBodyCenteredScene(const VECTOR3& tetColor = VECTOR3(0.5, 0.5, 0.5), 
                             const VECTOR3& outlineColor = VECTOR3(0,0,0))
  {
    cout << __FILE__ << " " << __FUNCTION__ << " " << __LINE__ << " : " << endl;
    cout << " NOT IMPLEMENTED. THAT'S WHY THE SCREEN IS BLANK." << endl;
  }
#endif
};

}

#endif
