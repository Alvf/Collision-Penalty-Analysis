#ifndef SHELL_SCENE_H
#define SHELL_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TRIANGLE_MESH.h"
#include "Timestepper/Shell/TIMESTEPPER.h"
#include "Hyperelastic/Shell/STVK.h"
#include "Hyperelastic/Shell/ARAP.h"
#include "Hyperelastic/Shell/BW_STRETCH.h"
#include "Hyperelastic/Shell/BW_SHEAR.h"
#include "Hyperelastic/Shell/BARAFF_WITKIN.h"
#include "Hyperelastic/Shell/BENDING_SPRING.h"
#include "Hyperelastic/Shell/DIHEDRAL.h"
#include "Collision/C_PLANES_EE.h"

namespace HOBAK {

class SHELL_SCENE : public SIMULATION_SCENE {
public:
  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Debugging shell scene" << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "shell";

    //_triangleMeshFilename = string("../data/objs/curtain_10x10.obj");
    _triangleMeshFilename = string("../data/objs/curtain_25x25.obj");
    //_triangleMeshFilename = string("../data/objs/curtain_50x50.obj");
    vector<VECTOR3> vertices;
    vector<VECTOR3I> triangles;
    bool success = TRIANGLE_MESH::readObjFile(_triangleMeshFilename, vertices, triangles);
    if (!success)
    {
      cout << " Failed to open file " << _triangleMeshFilename.c_str() << endl;
      return false;
    }
    _triangleMesh = new TRIANGLE_MESH(vertices, triangles);

    //_gravity = VECTOR3(0, 0, -1.0);
    _gravity = VECTOR3(-1.0, 0.0, 0.0);

    //_strechingEnergy = new SHELL::STVK(10.0, 1.0);
    _strechingEnergy = new SHELL::ARAP(10.0, 0.0);
    //_strechingEnergy = new SHELL::BW_STRETCH(1.0, 0.0);
    //_strechingEnergy = new SHELL::BW_SHEAR(1.0, 0.0);
    //_strechingEnergy = new SHELL::BARAFF_WITKIN(1.0, 0.0);
    //_bendingEnergy = new SHELL::BENDING_SPRING(20);
    //_bendingEnergy = new SHELL::BENDING_SPRING(10);
    _bendingEnergy = new SHELL::DIHEDRAL(0.5);
    // _bendingEnergy = new SHELL::DIHEDRAL(1);
    // _bendingEnergy = new SHELL::DIHEDRAL(10);
    //_bendingEnergy = new SHELL::DIHEDRAL(100);

    setCollisions(1000.0, 0.02);
    _shellSolver = new SHELL::TIMESTEPPER(*_triangleMesh, *_strechingEnergy, *_bendingEnergy, *_vfGeneral, *_eeGeneral);

    /*
    //Build the collision penalties
    vector<REAL> springArgs;
    springArgs.push_back(1000.0);
    springArgs.push_back(0.02);
    ENERGY_1D* normalSpring = new ENERGY_1D(springArgs);
    SIGNED_LEN_PLANES* lengthFunc = new SIGNED_LEN_PLANES();
    C_PLANES* cFunc = new C_PLANES();

    SIGNED_LEN_PLANES* lengthFuncEE = new SIGNED_LEN_PLANES();
    ENERGY_1D* normalSpringEE = new ENERGY_1D(springArgs);
    C_PLANES_EE* cFuncEE = new C_PLANES_EE();

    ENERGY_12D* vfGeneral = new ENERGY_12D(normalSpring, cFunc, lengthFunc);
    ENERGY_12D* eeGeneral = new ENERGY_12D(normalSpringEE, cFuncEE, lengthFuncEE);
    // build the time integrator
    _shellSolver = new SHELL::TIMESTEPPER(*_triangleMesh, *_strechingEnergy, *_bendingEnergy, *vfGeneral, *eeGeneral);
    //_shellSolver->setDt(1.0 / 60.0);
    */

    // cube on top
    VECTOR3 center0(0.0, 0.0, 0.95);
    //_kinematicShapes.push_back(new CUBE(center0, 1.0));
    //_shellSolver->attachKinematicSurfaceConstraints(_kinematicShapes[0]);
    //_shellSolver->addKinematicCollisionObject(_kinematicShapes[0]);
    
    VECTOR3 center1(-1.0, 0.02, -0);
    VECTOR3 center2(-5.0, 0.00, -0);
    //VECTOR3 center1(0.0, -1.0, -0.01);
    //_kinematicShapes.push_back(new SPHERE(center1, 0.5));
    _kinematicShapes.push_back(new SPHERE(center1, 0.25));
    _kinematicShapes.push_back(new CUBE(center2, 8));
    //_shellSolver->addKinematicCollisionObject(_kinematicShapes[1]);

   for(int i = 0; i < _kinematicShapes.size(); i++) 
    _shellSolver->addKinematicCollisionObject(_kinematicShapes[i]);

    //_eye    = VECTOR3(1.24999, 2.07096, 0.502227);
    //_lookAt = VECTOR3(0.777846, 1.20965, 0.314523);
    //_up     = VECTOR3(-0.0859995, -0.166908, 0.982213);

    _eye    = VECTOR3(2.03618, 0.26034, -1.50122);
    _lookAt = VECTOR3(1.26717, 0.000146806, -0.917327);
    _up     = VECTOR3(-0.190077, 0.965171, 0.179761);

    _worldCenter = _triangleMesh->getRestTranslation();
    //_worldCenter = VECTOR3(0, 0, 0);
    _pauseFrame = 37;

    // try scaling see if F follows
    //for (unsigned int x = 0; x < _triangleMesh->vertices().size(); x++)
    //  _triangleMesh->vertex(x)[2] *= 2.0;

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override {
    _shellSolver->externalForces().setZero();
    _shellSolver->addGravity(_gravity);
    _shellSolver->solve(verbose);

    //_kinematicShapes[1]->translation() += VECTOR3(0.0, 0.01, 0.0);
    //_kinematicShapes[1]->translation() += VECTOR3(0.0, 0.02, 0.0);
    _frameNumber++;
  };

// drawScene MUST have this guard, so it can be deactivated for regression testing
#ifndef GL_DISABLED
  // what to draw?
  virtual void drawScene()
  {
    drawAxes();
    drawTriangleMesh(*_triangleMesh, true);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);

    drawKinematicConstraints(_triangleMesh, _shellSolver);
    //drawPlaneConstraints(_triangleMesh, _shellSolver);

    drawVertex(*_triangleMesh, _arrowCounter);
    /*
    SIMULATION_SCENE::drawScene();

    if (_drawFeature)
    {
      glDisable(GL_DEPTH_TEST);
      glPointSize(10.0);
      drawVertexFacePairs(*_tetMesh, _arrowCounter);
      //drawVertexFacePairs(*_tetMesh);
    }
    */
  };
#endif
};

}

#endif
