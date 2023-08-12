#ifndef JSON_SHELL_SCENE_H
#define JSON_SHELL_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"

namespace HOBAK {

// Replaying simulation output from JSON file.
class JSON_SHELL_SCENE : public JSON_SCENE
{
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Replaying JSON output from shell simulateScene" << endl;
    cout << "=====================================================================" << endl;
  }

  // don't do anything for building the scene - FILE_IO should be doing
  // the heavy lifting here
  virtual bool buildScene() override
  {
    vector<VECTOR3> vertices;
    vector<VECTOR3I> triangles;
    bool success = TRIANGLE_MESH::readObjFile(_triangleMeshFilename, vertices, triangles);

    if (!success)
    {
      cout << " Failed to open file " << _triangleMeshFilename.c_str() << endl;
      return false;
    }

    if (_normalizedVertices)
      vertices = TRIANGLE_MESH::normalizeVertices(vertices);

    // apply any initial transforms the scene requested
    for (unsigned int x = 0; x < vertices.size(); x++)
      vertices[x] = _initialA * vertices[x] + _initialTranslation;

    _triangleMesh = new TRIANGLE_MESH(vertices, triangles);

    _drawFrame = 0;
    _drawFeature = true;
    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) {
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    
    _triangleMesh->setDisplacement(_positions[_drawFrame]);
    _triangleMesh->setCollisionPairs(_vertexFaceCollisions[_drawFrame],
                                     _edgeEdgeCollisions[_drawFrame]);

    if (verbose)
      cout << " Setting to frame " << _drawFrame << endl;
    _drawFrame++;
  };

  // jump to a specific frame
  void jumpToFrame(const int frame) {
    _drawFrame = frame;
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    
    _triangleMesh->setDisplacement(_positions[_drawFrame]);
    _triangleMesh->setCollisionPairs(_vertexFaceCollisions[_drawFrame],
                                     _edgeEdgeCollisions[_drawFrame]);

    cout << " Jumping to frame " << _drawFrame << endl;
  };

  // in case _drawFrame changed, you can update the positions
  void updatePositions()
  {
    if (_drawFrame >= (int)_positions.size()) return;
    if (_drawFrame < 0) return;
    _triangleMesh->setDisplacement(_positions[_drawFrame]);
    _triangleMesh->setCollisionPairs(_vertexFaceCollisions[_drawFrame],
                                     _edgeEdgeCollisions[_drawFrame]);
  }

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene() override
  {
    glEnable(GL_DEPTH_TEST);
    drawTriangleMesh(*_triangleMesh, true);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
   
    glEnable(GL_DEPTH_TEST);
    // TODO: implement for triangle meshes
    //if (_drawFeature)
    //  drawVertexFacePairs(*_triangleMesh, _arrowCounter);
  };

  // subtract off the translation and rotation from the main body
  // and draw that
  void drawBodyCenteredScene(const VECTOR3& tetColor = VECTOR3(0.5, 0.5, 0.5), 
                             const VECTOR3& outlineColor = VECTOR3(0,0,0))
  {
    glEnable(GL_DEPTH_TEST);

    // let's orient ourselves
    drawAxes();

    glPushMatrix();

    // undo the global rotation
    const MATRIX3 RT = _triangleMesh->getRotation().transpose();
    const Eigen::AngleAxis<GLfloat> rotation{ RT.cast<GLfloat>() };
    glRotatef((180.0 / M_PI) * rotation.angle(), rotation.axis().x(), rotation.axis().y(), rotation.axis().z());
    
    // undo the global translation
    const VECTOR3 translation = _triangleMesh->getTranslation();
    glTranslatef(-translation[0], -translation[1], -translation[2]);

    // now draw
    drawTriangleMesh(*_triangleMesh, true);

    // TODO: implement for triangle meshes
    //if (_drawFeature)
    //  drawVertexFacePairs(*_triangleMesh, _arrowCounter);

    glPopMatrix();
  }
#endif
};

}

#endif
