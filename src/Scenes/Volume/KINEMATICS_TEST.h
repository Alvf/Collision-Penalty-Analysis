#ifndef KINEMATICS_TEST_H
#define KINEMATICS_TEST_H

#include "Scenes/SIMULATION_SCENE.h"

namespace HOBAK {

class KINEMATICS_TEST : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Smashing a cube between two kinematic colliders, to see if " << endl;
    cout << " kinematic constraints are working correctly" << endl;
    cout << endl;
    cout << " Left and right arrow keys move the left cube back and forth" << endl;
    cout << " You can see the cube react accordingly. " << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "kinematics_test";
    // read in the tet mesh file
    _tetMeshFilename = string("../data/cube_6.tobj");
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
    _hyperelastic = new VOLUME::ARAP(1.0, 1.0);

    setCollisions();

    // build the time integrator
    _volumeSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);

    // add the kinematics
    const VECTOR3 center0(0.56,0.48, -0.25);
    const VECTOR3 center1(0.56,0.48,1.23);
    _kinematicShapes.push_back(new CUBE(center0, 1.0));
    _kinematicShapes.push_back(new CUBE(center1, 1.0));

    // attach it to the rest of the scene
    _volumeSolver->attachKinematicSurfaceConstraints(_kinematicShapes[1]);
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes[0]);

    _eye    = VECTOR3(-1.61233, 0.617816, -0.169647);
    _lookAt = VECTOR3(-0.630224, 0.573968, 0.0135059);
    _up     = VECTOR3(0.0418256, 0.999014, 0.0148915);
    return true;
  }
 
  // cycle the cube according to a scripted motion
  void cubeMotion()
  {
    // TODO: this won't work because the velocity of each kinematic object
    // is not currently stored. The relative velocity when doing collision detection
    // in findNewSurfaceConstraints will then show up as zero, and it will not
    // form a new plane constraint.
    //
    //const int whichCube = 0;
    
    // currently this one works because by the time the cube pushes into the
    // other plane, it has a velocity, and the relative velocity show up
    // as non-zero
    const int whichCube = 1;

    if (_frameNumber > 1 && _frameNumber < 7)
    {
      _kinematicShapes[whichCube]->translation()[2] -= 0.01;
      return;
    }
    
    if (_frameNumber >= 7 && _frameNumber < 13)
    {
      _kinematicShapes[whichCube]->translation()[2] += 0.01;
      return;
    }
    
    if (_frameNumber >= 17 && _frameNumber < 30)
    {
      _kinematicShapes[whichCube]->translation()[2] -= 0.01;
      return;
    }
    
    if (_frameNumber >= 30 && _frameNumber < 45)
    {
      _kinematicShapes[whichCube]->translation()[2] += 0.01;
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
      _kinematicShapes[1]->translation()[2] -= 0.01;
      _leftArrow = false;
    }
    if (_rightArrow)
    {
      _kinematicShapes[1]->translation()[2] += 0.01;
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
