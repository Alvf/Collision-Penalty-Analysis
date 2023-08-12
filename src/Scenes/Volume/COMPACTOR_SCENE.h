#ifndef COMPACTOR_SCENE_H
#define COMPACTOR_SCENE_H

#include "Scenes/SIMULATION_SCENE.h"
#include "Timestepper/Volume/BDF_1.h"
#include "Timestepper/Volume/BDF_2.h"
#include "Collision/C_PLANES_EE.h"

namespace HOBAK {

class COMPACTOR_SCENE : public SIMULATION_SCENE {
public:

  virtual void printSceneDescription() override
  {
    cout << "=====================================================================" << endl;
    cout << " Compacting a bunny for the SCA paper revision" << endl;
    cout << "=====================================================================" << endl;
  }

  virtual bool buildScene() override
  {
    // this will determine the MOV and JSON filenames
    _sceneName = "bunny_compactor";

    // read in the tet mesh file
    //_tetMeshFilename = string("../data/scorpion_0_125.tobj");
    // _tetMeshFilename = string("../data/bunny/bunny_5.tobj");
    string cubeMesh = string("../data/single_tet.tobj");
    //_tetMeshFilename = string("../data/Brachiosaurus_15.tobj");
    //_tetMeshFilename = string("../data/dragon.tobj");
    vector<VECTOR3> vertices;
    vector<VECTOR4I> tets;
    //the tobjs are REALLY WHACK; doing normalizeVerts after
    //these translates means you have to fiddle here...
    bool success = TET_MESH::readTobjFile(cubeMesh, vertices, tets, false, VECTOR3(0.0,0.0,0.0));
    // success = success && TET_MESH::readTobjFile(cubeMesh, vertices, tets, true, VECTOR3(,,));
    if (!success)
    {
      cout << " Failed to open file " << _tetMeshFilename.c_str() << endl;
      return false;
    }
    vertices = TET_MESH::normalizeVertices(vertices);
    // cout << "skipped normalizeVertices" << endl;
    _normalizedVertices = true;

    using namespace Eigen;

    _gravity = VECTOR3(0, -1.0, 0);

    REAL E = 5.0;
    REAL nu = 0.45; // lambda \approx 10

    REAL mu     = VOLUME::HYPERELASTIC::computeMu(E, nu);
    REAL lambda = VOLUME::HYPERELASTIC::computeLambda(E, nu);
    cout << " mu:     " << mu << endl;
    cout << " lambda: " << lambda << endl;

    // build the tet mesh object
    _tetMesh = new TET_MESH_FASTER(vertices, tets);
    _hyperelastic = new VOLUME::SNH(mu, lambda);
    //_hyperelastic = new VOLUME::ARAP(mu * 4, lambda);

    const vector<REAL>& areas = _tetMesh->surfaceTriangleAreas();
    REAL smallest = areas[0];
    REAL largest = areas[0];
    for (unsigned int x = 1; x < areas.size(); x++)
    {
      if (areas[x] > largest) largest  = areas[x];
      if (areas[x] < largest) smallest = areas[x];
    }
    cout << "Largest triangle area:  "  << largest << endl;
    cout << "Smallest triangle area: "  << smallest << endl;

    setCollisions(1200.0, 0.02);
    _volumeSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic, *_vfGeneral, *_eeGeneral);

    /*
    vector<REAL> springArgs;
    springArgs.push_back(1200.0);
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
    _volumeSolver = new VOLUME::BACKWARD_EULER_POSITION(*_tetMesh, *_hyperelastic, *vfGeneral, *eeGeneral);
    */

    //_volumeSolver = new VOLUME::BDF_2(*_tetMesh, *_hyperelastic);
    //
    // FLT_EPSILON in VOLUME::TIMESTEPPER::findNewSurfaceConstraints does not play nice here;
    // it seems to like it when it's set to zero. Needs further investigation.
    // Plus the conditioning of the matrix seems much worse than position-based
    //_volumeSolver = new VOLUME::BACKWARD_EULER_VELOCITY(*_tetMesh, *_hyperelastic);
    _volumeSolver->setDt(1.0 / 40.0);

    _allCenter= VECTOR3(1.0,1.0,1.0);

    VECTOR3 center(0.0, -2, 0.0);
    _centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(_centers.back() + _allCenter, 2.0));
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());
    
    center = VECTOR3(0, 2.0, 0);
    _centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(_centers.back() + _allCenter, 2.0));
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(2.0, 0.0, 0.0);
    _centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(_centers.back() + _allCenter, 2.0));
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(-2.0, 0.0, 0.0);
    _centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(_centers.back() + _allCenter, 2.0));
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(0.0, 0.0, 2.0);
    _centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(_centers.back() + _allCenter, 2.0));
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());

    center = VECTOR3(0.0, 0.0, -2.00);
    _centers.push_back(center);
    _kinematicShapes.push_back(new CUBE(_centers.back() + _allCenter, 2.0));
    _volumeSolver->addKinematicCollisionObject(_kinematicShapes.back());
    // collision constants
    const REAL collisionMu = 1000.0;
    _volumeSolver->collisionStiffness() = collisionMu;
    _volumeSolver->collisionDampingBeta() = 0.01;
    
    _volumeSolver->vertexFaceSelfCollisionsOn() = true;
    _volumeSolver->edgeEdgeSelfCollisionsOn() = true;

    _eye    = VECTOR3(6.18475, 1.28356, 2.91415);
    _lookAt = VECTOR3(5.26567, 1.23173, 2.52334);
    _up     = VECTOR3(-0.211152, 0.901845, 0.37695);

    _worldCenter = VECTOR3( 0.497, 0.785889, 0.452556);
    //_pauseFrame = 800;
    _pauseFrame = 400;

    return true;
  }
  virtual void stepSimulation(const bool verbose = true) {
    _volumeSolver->externalForces().setZero();
    _volumeSolver->addGravity(_gravity);
    _volumeSolver->solve(verbose);

    if (_writeToFile) {
      if (!writeFrameToFile()) {
        // Do something here that isn't too noisy.
      }
    }

    //compression
    if (_frameNumber <= _maxCompressFrame){
      for(unsigned int i = 0; i < _kinematicShapes.size(); i++){
        _kinematicShapes[i]->translation() = _centers[i]*(1.0 - _compressionFactor*_frameNumber/_maxCompressFrame) + _allCenter; 
      }
    }

    //coming back
    if (_frameNumber >= _startOpenFrame && _frameNumber <= _maxCompressFrame){
      for(unsigned int i = 0; i < _kinematicShapes.size(); i++){
        _kinematicShapes[i]->translation() = _centers[i]*(1.0 - _compressionFactor*_frameNumber/_maxCompressFrame) + _allCenter; 
      }
    }


    _frameNumber++;
  };

  virtual void drawScene() override{
    glEnable(GL_DEPTH_TEST);
    drawSurfaceTriangles(*_tetMesh, true);
    for (unsigned int x = 0; x < _kinematicShapes.size(); x++)
      drawKinematicShape(*_kinematicShapes[x]);
    drawAxes();
  }
private:
  vector<VECTOR3> _centers;
  VECTOR3 _allCenter;
  REAL _compressionFactor = 1.0/2.76;
  int _maxCompressFrame = 200;
  int _startOpenFrame = 300;
};

}

#endif
