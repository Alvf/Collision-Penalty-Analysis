#ifndef NEWMARK_STRETCH_H
#define NEWMARK_STRETCH_H

#include "Scenes/SIMULATION_SCENE.h"

namespace HOBAK {

class NEWMARK_STRETCH : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Stretch out a cube with a Newmark solver, see it wiggle" << endl;
    cout << endl;
    cout << " You can move the left cube with left and right arrow keys " << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "newmark_stretch";

    // read in the tet mesh file
    _tetMeshFilename = string("../data/cube_10.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    bool success = TET_MESH::readTobjFile(_tetMeshFilename, vertices, tets);

    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }

    // build the tet mesh object
    _tetMesh = new TET_MESH(vertices, tets);
    _hyperelastic = new VOLUME::SNH(1.0, 10.0);
    
    // setup collision objects
    setCollisions(1.0, 0.01);

    // build the time integrator
    _volumeSolver = new VOLUME::NEWMARK(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);

    // add the kinematics
    const VECTOR3 center0(0.56,0.48, -0.24);
    const VECTOR3 center1(0.56,0.48,1.23);
    _kinematicShapes.push_back(new CUBE(center0, 1.0));
    _kinematicShapes.push_back(new CUBE(center1, 1.0));

    // attach it to the rest of the scene
    _volumeSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    _volumeSolver->attachKinematicSurfaceConstraints(_kinematicShapes[1]);

    _eye    = VECTOR3(-1.73357, 0.953765, 0.471698);
    _lookAt = VECTOR3(-0.749345, 0.790722, 0.402909);
    _up     = VECTOR3(0.160779, 0.986286, -0.0372856);

    _pauseFrame = 75;
    return true;
  }
 
  // cycle the cube according to a scripted motion
  void cubeMotion()
  {
    if (_frameNumber > 1 && _frameNumber < 25)
    {
      _kinematicShapes[0]->translation()[2] -= 0.01;
      return;
    }
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override
  {
    // scripted motion
    cubeMotion();

    _volumeSolver->externalForces().setZero();
    _volumeSolver->addGravity(_gravity);
    _volumeSolver->solve(verbose);

    if (_leftArrow)
    {
      _kinematicShapes[0]->translation()[2] -= 0.01;
      _leftArrow = false;
    }
    if (_rightArrow)
    {
      _kinematicShapes[0]->translation()[2] += 0.01;
      _rightArrow = false;
    }

    _frameNumber++;
  }

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene() override
  {
    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    glDisable(GL_DEPTH_TEST);
    glPointSize(10.0);
    drawVertices(*_tetMesh, _volumeSolver->constrainedNodes());
    drawPlaneConstraints(_tetMesh, _volumeSolver);
  };
#endif

};

}

#endif
