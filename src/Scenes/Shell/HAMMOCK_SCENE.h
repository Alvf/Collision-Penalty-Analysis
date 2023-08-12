#ifndef HAMMOCK_SCENE_H
#define HAMMOCK_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Geometry/TRIANGLE_MESH.h"
#include "Geometry/TRIANGLE_MESH_FASTER.h"
#include "Timestepper/Shell/TIMESTEPPER.h"
#include "Hyperelastic/Shell/ARAP.h"
#include "Hyperelastic/Shell/DIHEDRAL.h"
#include "Collision/ENERGY_12D.h"
#include "Collision/C_PLANES_EE.h"

namespace HOBAK {

std::vector<KINEMATIC_SHAPE*> cubes;
vector<int> dropFrames;
REAL dt = 1/40.0;
REAL dropInterval = 3.25;
int frameInterval = dropInterval/dt;
SHELL::TIMESTEPPER* clothSolver;

class HAMMOCK_SCENE : public SIMULATION_SCENE {
public:
  HAMMOCK_SCENE(){
  }

  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Debugging hammock scene" << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    _sceneName = "Hammock_Scene";

    std::string hammockFile = string("../data/objs/HammockTilt.obj");
    std::string hammockDrop1 = string("../data/objs/HammockDrop1.obj");
    std::string hammockDrop2 = string("../data/objs/HammockDrop2.obj");
    std::string hammockDrop3 = string("../data/objs/HammockDrop3.obj");

    vector<VECTOR3> vertices;
    vector<VECTOR3I> triangles;
    bool success = TRIANGLE_MESH::readObjFile(hammockFile, vertices, triangles,false);

    vector<std::string> dropSubjects;
    dropSubjects.push_back(hammockDrop1);
    dropSubjects.push_back(hammockDrop2);
    dropSubjects.push_back(hammockDrop3);

    // success = success && TRIANGLE_MESH::readObjFile(mesh2, vertices, triangles,true);
    VECTOR3 translator = VECTOR3(3,3,0);
    for(unsigned int i = 0; i < dropSubjects.size(); i++){ //loop limit here determines how many sheets get dropped after the initial one
      cubes.push_back(new CUBE(translator,2.0));
      success = success && TRIANGLE_MESH::readObjFile(dropSubjects[i], vertices, triangles,true, translator);
      if (!success)
      {
        cout << " Failed to open file " << dropSubjects[i] << endl;
        return false;
      }
      translator[0] += 2;
    }
    dropFrames.push_back(20);
    for(unsigned int i = 0; i < cubes.size() - 1; i ++) dropFrames.push_back(frameInterval*(1+i));
    _triangleMesh = new TRIANGLE_MESH_FASTER(vertices, triangles);

    _gravity = VECTOR3(0, -1, 0.0);

    // building the energies
    _strechingEnergy = new SHELL::ARAP(10.0, 0.0);
    _bendingEnergy = new SHELL::DIHEDRAL(0.01);
    // swap these out to pick between unsigned and signed collisions
    // setCollisions(1000.0, 0.02);
    setUnsignedCollisions(1000.0, 0.02);
    // build the time integrator
    clothSolver = new SHELL::TIMESTEPPER(*_triangleMesh, *_strechingEnergy, *_bendingEnergy, *_vfGeneral, *_eeGeneral);
    clothSolver->setDt(dt);

    //pinning cubes (bottom)
    VECTOR3 c10(0.45, -0.325, -0.25);
    VECTOR3 c11(0.45, -0.325, 0.25);
    VECTOR3 c12(-0.45, -0.03, -0.25);
    VECTOR3 c13(-0.45, -0.03, 0.25);
    _kinematicShapes.push_back(new CUBE(c10, 0.05));
    _kinematicShapes.push_back(new CUBE(c11, 0.05));
    _kinematicShapes.push_back(new CUBE(c12, 0.05));
    _kinematicShapes.push_back(new CUBE(c13, 0.05));
    for(int i = 0; i < _kinematicShapes.size(); i++){
      clothSolver->attachKinematicSurfaceConstraints(_kinematicShapes[i]);
    }
    for (int i = 0; i < cubes.size(); i++) { // set the pinning cubes
      clothSolver->attachKinematicSurfaceConstraints(cubes[i]);
    }
    clothSolver->collisionDampingBeta() = 0;

    _eye    = VECTOR3(0.72963, 0.456088, 0.813303); //closer to the action
    _lookAt = VECTOR3(0.199469, -0.133773, 0.204212);
    _up     = VECTOR3(-0.365711, 0.807177, -0.463374);

    _worldCenter = _triangleMesh->getRestTranslation();
    _pauseFrame = 500; //long enough for 4 dropped layers 

    return true;
  }

  // simulation loop
  virtual void stepSimulation(const bool verbose = true) override {
    for(int i = 0; i < dropFrames.size(); i++){
      if(_frameNumber == dropFrames[i]){
        cout << "teleporting drop cube on frame" << dropFrames[i] << endl;
        cubes[i]->translation() = VECTOR3(0,0,0);
      }
    }

    clothSolver->externalForces().setZero();
    clothSolver->addGravity(_gravity);
    clothSolver->solve(verbose);
    
    for(int i = 0; i < dropFrames.size(); i++){
      if(_frameNumber == dropFrames[i]){
        clothSolver->releaseKinematicSurfaceConstraints(cubes[i]);
      }
    }

    _frameNumber++;
  };

  virtual void drawScene()
  {
    drawAxes();
    drawTriangleMesh(*_triangleMesh, true);

    glEnable(GL_DEPTH_TEST);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);


    drawVertex(*_triangleMesh, _arrowCounter);
  };

protected:
  TRIANGLE_MESH* _triangleMesh;
  SHELL::STRETCHING* _strechingEnergy;
  SHELL::BENDING* _bendingEnergy;

  string _triangleMeshFilename;
};

}

#endif
